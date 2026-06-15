from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


@dataclass(frozen=True)
class ProgramMatch:
    reference_program: str
    candidate_program: str
    spectra_cosine: float
    usage_pearson: float
    top_gene_jaccard: float


@dataclass(frozen=True)
class KConsistency:
    k: int
    matches: tuple[ProgramMatch, ...]
    mean_spectra_cosine: float
    min_spectra_cosine: float
    mean_usage_pearson: float
    min_usage_pearson: float
    mean_top_gene_jaccard: float
    overall_consistency: float


@dataclass(frozen=True)
class BenchmarkGateResult:
    reference_factorize_seconds: float
    candidate_factorize_seconds: float
    speedup: float
    time_saved_fraction: float
    min_required_time_saved_fraction: float
    min_required_consistency: float
    consistency: tuple[KConsistency, ...]
    passed: bool

    def to_json(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2), encoding="utf-8")


def read_cnmf_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def _cosine_similarity_rows(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a_norm = np.linalg.norm(a, axis=1, keepdims=True)
    b_norm = np.linalg.norm(b, axis=1, keepdims=True)
    a_safe = np.divide(a, a_norm, out=np.zeros_like(a, dtype=float), where=a_norm != 0)
    b_safe = np.divide(b, b_norm, out=np.zeros_like(b, dtype=float), where=b_norm != 0)
    return a_safe @ b_safe.T


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a - a.mean()
    b = b - b.mean()
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom == 0:
        return 0.0
    return float(np.dot(a, b) / denom)


def _top_gene_jaccard(ref: pd.Series, cand: pd.Series, top_n: int = 50) -> float:
    top_ref = set(ref.sort_values(ascending=False).head(top_n).index)
    top_cand = set(cand.sort_values(ascending=False).head(top_n).index)
    union = top_ref | top_cand
    if not union:
        return 1.0
    return len(top_ref & top_cand) / len(union)


def compare_cnmf_k(
    reference_spectra: Path,
    candidate_spectra: Path,
    reference_usage: Path,
    candidate_usage: Path,
    k: int,
    top_gene_n: int = 50,
) -> KConsistency:
    ref_spectra = read_cnmf_matrix(reference_spectra)
    cand_spectra = read_cnmf_matrix(candidate_spectra)
    ref_usage = read_cnmf_matrix(reference_usage)
    cand_usage = read_cnmf_matrix(candidate_usage)

    genes = ref_spectra.columns.intersection(cand_spectra.columns)
    if len(genes) == 0:
        raise ValueError(f"no overlapping genes for k={k}")
    obs = ref_usage.index.intersection(cand_usage.index)
    if len(obs) == 0:
        raise ValueError(f"no overlapping observations for k={k}")

    ref_s = ref_spectra.loc[:, genes].astype(float)
    cand_s = cand_spectra.loc[:, genes].astype(float)
    cosine = _cosine_similarity_rows(ref_s.to_numpy(), cand_s.to_numpy())
    ref_idx, cand_idx = linear_sum_assignment(-cosine)

    matches: list[ProgramMatch] = []
    for r_i, c_i in zip(ref_idx, cand_idx):
        ref_program = str(ref_s.index[r_i])
        cand_program = str(cand_s.index[c_i])
        ref_usage_col = ref_program if ref_program in ref_usage.columns else str(int(float(ref_program)))
        cand_usage_col = cand_program if cand_program in cand_usage.columns else str(int(float(cand_program)))
        usage_corr = _pearson(
            ref_usage.loc[obs, ref_usage_col].to_numpy(),
            cand_usage.loc[obs, cand_usage_col].to_numpy(),
        )
        jaccard = _top_gene_jaccard(ref_s.iloc[r_i], cand_s.iloc[c_i], top_n=top_gene_n)
        matches.append(
            ProgramMatch(
                reference_program=ref_program,
                candidate_program=cand_program,
                spectra_cosine=float(cosine[r_i, c_i]),
                usage_pearson=usage_corr,
                top_gene_jaccard=float(jaccard),
            )
        )

    spectra_values = np.array([m.spectra_cosine for m in matches], dtype=float)
    usage_values = np.array([m.usage_pearson for m in matches], dtype=float)
    jaccard_values = np.array([m.top_gene_jaccard for m in matches], dtype=float)
    overall = min(float(spectra_values.mean()), float(usage_values.mean()))
    return KConsistency(
        k=k,
        matches=tuple(matches),
        mean_spectra_cosine=float(spectra_values.mean()),
        min_spectra_cosine=float(spectra_values.min()),
        mean_usage_pearson=float(usage_values.mean()),
        min_usage_pearson=float(usage_values.min()),
        mean_top_gene_jaccard=float(jaccard_values.mean()),
        overall_consistency=overall,
    )


def cnmf_output_paths(output_dir: Path, run_name: str, k: int) -> tuple[Path, Path]:
    spectra = output_dir / f"{run_name}.spectra.k_{k}.dt_0_5.consensus.txt"
    usage = output_dir / f"{run_name}.usages.k_{k}.dt_0_5.consensus.txt"
    return spectra, usage


def compare_cnmf_dirs(
    *,
    reference_dir: Path,
    candidate_dir: Path,
    reference_run_name: str,
    candidate_run_name: str,
    k_values: tuple[int, ...],
    output_json: Path,
) -> dict:
    consistency = []
    for k in k_values:
        ref_spectra, ref_usage = cnmf_output_paths(reference_dir, reference_run_name, k)
        cand_spectra, cand_usage = cnmf_output_paths(candidate_dir, candidate_run_name, k)
        consistency.append(
            compare_cnmf_k(
                reference_spectra=ref_spectra,
                candidate_spectra=cand_spectra,
                reference_usage=ref_usage,
                candidate_usage=cand_usage,
                k=k,
            )
        )
    mean_spectra = float(np.mean([item.mean_spectra_cosine for item in consistency]))
    mean_usage = float(np.mean([item.mean_usage_pearson for item in consistency]))
    min_k_overall = float(min(item.overall_consistency for item in consistency))
    result = {
        "reference_dir": str(reference_dir),
        "candidate_dir": str(candidate_dir),
        "reference_run_name": reference_run_name,
        "candidate_run_name": candidate_run_name,
        "k_values": list(k_values),
        "mean_spectra_cosine": mean_spectra,
        "min_spectra_cosine": float(min(item.min_spectra_cosine for item in consistency)),
        "mean_usage_pearson": mean_usage,
        "min_usage_pearson": float(min(item.min_usage_pearson for item in consistency)),
        "passes_95pct_gate": bool(min(mean_spectra, mean_usage) >= 0.95),
        "min_k_overall_consistency": min_k_overall,
        "passes_all_k_95pct_gate": bool(min_k_overall >= 0.95),
        "comparisons": [asdict(item) for item in consistency],
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def load_factorize_times(summary_json: Path) -> tuple[float, float]:
    summary = json.loads(summary_json.read_text(encoding="utf-8"))
    serial = float(summary["comparison"]["serial_factorize_seconds"])
    parallel = float(summary["comparison"]["parallel_factorize_seconds"])
    return serial, parallel


def evaluate_benchmark_gate(
    summary_json: Path,
    reference_dir: Path,
    candidate_dir: Path,
    run_name: str,
    k_values: tuple[int, ...],
    min_time_saved_fraction: float = 0.50,
    min_consistency: float = 0.95,
) -> BenchmarkGateResult:
    reference_time, candidate_time = load_factorize_times(summary_json)
    consistency: list[KConsistency] = []
    for k in k_values:
        ref_spectra, ref_usage = cnmf_output_paths(reference_dir, run_name, k)
        cand_spectra, cand_usage = cnmf_output_paths(candidate_dir, run_name, k)
        consistency.append(
            compare_cnmf_k(
                reference_spectra=ref_spectra,
                candidate_spectra=cand_spectra,
                reference_usage=ref_usage,
                candidate_usage=cand_usage,
                k=k,
            )
        )

    speedup = reference_time / candidate_time
    saved_fraction = (reference_time - candidate_time) / reference_time
    passed = (
        saved_fraction >= min_time_saved_fraction
        and all(kc.overall_consistency >= min_consistency for kc in consistency)
    )
    return BenchmarkGateResult(
        reference_factorize_seconds=reference_time,
        candidate_factorize_seconds=candidate_time,
        speedup=float(speedup),
        time_saved_fraction=float(saved_fraction),
        min_required_time_saved_fraction=min_time_saved_fraction,
        min_required_consistency=min_consistency,
        consistency=tuple(consistency),
        passed=passed,
    )
