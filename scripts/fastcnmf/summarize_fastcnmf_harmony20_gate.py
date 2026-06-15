#!/usr/bin/env python
"""Summarize FastCNMF Harmony2 end-to-end speed and consistency gate."""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


ROOT = Path(__file__).resolve().parents[2]
TMP = ROOT / "tmp/fastcnmf_harmony2"
BASE = ROOT / "tmp/cnmf_spatial_harmony_benchmark"
OUT_JSON = TMP / "fastcnmf_harmony20_gate.json"
OUT_MD = TMP / "fastcnmf_harmony20_gate.md"


def parse_elapsed(value: str) -> float:
    parts = value.strip().split(":")
    if len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    return float(value)


def parse_time(path: Path) -> dict:
    data = {"path": str(path)}
    fields = {
        "Elapsed (wall clock) time (h:mm:ss or m:ss):": ("elapsed_seconds", parse_elapsed),
        "Maximum resident set size (kbytes):": ("max_rss_mb", lambda x: round(int(x) / 1024, 1)),
        "Percent of CPU this job got:": ("cpu_percent", lambda x: float(x.rstrip("%"))),
        "Exit status:": ("exit_status", int),
    }
    for line in path.read_text(errors="replace").splitlines():
        line = line.strip()
        for prefix, (key, fn) in fields.items():
            if line.startswith(prefix):
                data[key] = fn(line.removeprefix(prefix).strip())
                break
    return data


def parallel_group_elapsed(path: Path) -> int:
    text = path.read_text()
    return int(re.search(r"elapsed_seconds=(\d+)", text).group(1))


def cosine_matrix(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    a_norm = np.linalg.norm(a, axis=1, keepdims=True)
    b_norm = np.linalg.norm(b, axis=1, keepdims=True)
    return (a @ b.T) / np.maximum(a_norm, 1e-12) / np.maximum(b_norm.T, 1e-12)


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def compare_k(k: int) -> dict:
    base_dir = BASE / "parallel/gbm_lowres_harmony_cnmf"
    fast_dir = TMP / "parallel/gbm_lowres_harmony20_fastcnmf"
    base_spectra = pd.read_csv(base_dir / f"gbm_lowres_harmony_cnmf.spectra.k_{k}.dt_0_5.consensus.txt", sep="\t", index_col=0)
    fast_spectra = pd.read_csv(fast_dir / f"gbm_lowres_harmony20_fastcnmf.spectra.k_{k}.dt_0_5.consensus.txt", sep="\t", index_col=0)
    common_genes = base_spectra.columns.intersection(fast_spectra.columns)
    base_s = base_spectra.loc[:, common_genes].to_numpy()
    fast_s = fast_spectra.loc[:, common_genes].to_numpy()
    sim = cosine_matrix(base_s, fast_s)
    rows, cols = linear_sum_assignment(-sim)

    base_usage = pd.read_csv(base_dir / f"gbm_lowres_harmony_cnmf.usages.k_{k}.dt_0_5.consensus.txt", sep="\t", index_col=0)
    fast_usage = pd.read_csv(fast_dir / f"gbm_lowres_harmony20_fastcnmf.usages.k_{k}.dt_0_5.consensus.txt", sep="\t", index_col=0)
    common_obs = base_usage.index.intersection(fast_usage.index)
    usage_corrs = []
    spectra_cos = []
    for r, c in zip(rows, cols):
        spectra_cos.append(float(sim[r, c]))
        usage_corrs.append(pearson(base_usage.loc[common_obs].iloc[:, r].to_numpy(), fast_usage.loc[common_obs].iloc[:, c].to_numpy()))
    return {
        "k": k,
        "common_genes": int(len(common_genes)),
        "common_obs": int(len(common_obs)),
        "program_match": [{"baseline": int(r + 1), "fast": int(c + 1), "spectra_cosine": float(sim[r, c])} for r, c in zip(rows, cols)],
        "mean_spectra_cosine": float(np.nanmean(spectra_cos)),
        "min_spectra_cosine": float(np.nanmin(spectra_cos)),
        "mean_usage_pearson": float(np.nanmean(usage_corrs)),
        "min_usage_pearson": float(np.nanmin(usage_corrs)),
    }


def main() -> None:
    baseline = json.loads((BASE / "benchmark_summary.json").read_text())
    original_serial = baseline["prepare_input"]["elapsed_seconds"] + sum(
        baseline["runs"]["serial"]["steps"][step]["elapsed_seconds"]
        for step in ["prepare", "factorize", "combine", "consensus", "k_selection_plot"]
    )

    fast_steps = {
        "core_prepare": parse_time(TMP / "prepare_core.time.log"),
        "harmony2_adapter": parse_time(TMP / "harmony20.time.log"),
        "materialize_input": parse_time(TMP / "materialize_input.time.log"),
        "cnmf_prepare": parse_time(TMP / "logs/parallel/prepare.time.log"),
        "factorize_parallel": {"elapsed_seconds": parallel_group_elapsed(TMP / "logs/parallel/factorize_parallel.cmd.log")},
        "combine": parse_time(TMP / "logs/parallel/combine.time.log"),
        "consensus": parse_time(TMP / "logs/parallel/consensus.time.log"),
        "k_selection_plot": parse_time(TMP / "logs/parallel/k_selection_plot.time.log"),
    }
    fast_total = sum(step["elapsed_seconds"] for step in fast_steps.values())
    comparisons = [compare_k(6), compare_k(8)]
    min_spectra = min(c["min_spectra_cosine"] for c in comparisons)
    min_usage = min(c["min_usage_pearson"] for c in comparisons)
    mean_spectra = float(np.mean([c["mean_spectra_cosine"] for c in comparisons]))
    mean_usage = float(np.mean([c["mean_usage_pearson"] for c in comparisons]))

    summary = {
        "baseline_original_cnmf_serial_end_to_end_seconds": original_serial,
        "fastcnmf_harmony20_end_to_end_seconds": fast_total,
        "overall_speedup": original_serial / fast_total,
        "passes_2x_speed_gate": original_serial / fast_total >= 2.0,
        "consistency": {
            "comparisons": comparisons,
            "mean_spectra_cosine": mean_spectra,
            "min_spectra_cosine": min_spectra,
            "mean_usage_pearson": mean_usage,
            "min_usage_pearson": min_usage,
            "passes_95pct_gate": min(mean_spectra, mean_usage) >= 0.95,
        },
        "fast_steps": fast_steps,
        "notes": [
            "FastCNMF uses harmonypy 2.0 C++ Harmony plus reconstructed Phi_moe and cNMF-compatible fixed lamb=1 MOE correction.",
            "The original cNMF baseline used cnmf 1.7.1 + harmonypy 0.2.0 preprocessing; harmonypy 0.2 defaults match cNMF's fixed-lambda MOE contract.",
        ],
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    md = f"""# FastCNMF Harmony2 Gate

## Speed

- Original cNMF serial end-to-end: `{original_serial:.2f} s`
- FastCNMF Harmony2 end-to-end: `{fast_total:.2f} s`
- Overall speedup: `{summary['overall_speedup']:.2f}x`
- Passes >2x overall speed gate: `{summary['passes_2x_speed_gate']}`

## Consistency

- Mean spectra cosine: `{mean_spectra:.4f}`
- Minimum spectra cosine: `{min_spectra:.4f}`
- Mean usage Pearson: `{mean_usage:.4f}`
- Minimum usage Pearson: `{min_usage:.4f}`
- Passes >95% mean consistency gate: `{summary['consistency']['passes_95pct_gate']}`

## Step Times

| step | seconds | max_rss_mb |
| --- | ---: | ---: |
"""
    for name, data in fast_steps.items():
        md += f"| {name} | {data.get('elapsed_seconds', math.nan)} | {data.get('max_rss_mb', '')} |\n"
    md += "\n## Per-k Consistency\n\n"
    for comp in comparisons:
        md += (
            f"- k={comp['k']}: mean spectra cosine `{comp['mean_spectra_cosine']:.4f}`, "
            f"mean usage Pearson `{comp['mean_usage_pearson']:.4f}`\n"
        )
    md += "\n## Notes\n\n"
    for note in summary["notes"]:
        md += f"- {note}\n"
    OUT_MD.write_text(md, encoding="utf-8")
    print(json.dumps({
        "speedup": summary["overall_speedup"],
        "passes_speed": summary["passes_2x_speed_gate"],
        "mean_spectra_cosine": mean_spectra,
        "mean_usage_pearson": mean_usage,
        "passes_consistency": summary["consistency"]["passes_95pct_gate"],
    }, indent=2))


if __name__ == "__main__":
    main()
