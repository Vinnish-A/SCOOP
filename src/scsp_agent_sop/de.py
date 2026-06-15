from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
from scipy import stats


@dataclass(frozen=True)
class PseudobulkDEResult:
    table: pd.DataFrame
    design: pd.DataFrame
    manifest: dict


def read_pseudobulk_dir(path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = Path(path)
    counts = pd.read_csv(path / "counts.tsv", sep="\t")
    meta = pd.read_csv(path / "metadata.tsv", sep="\t")
    if "pseudobulk_id" not in counts.columns or "pseudobulk_id" not in meta.columns:
        raise ValueError("counts.tsv and metadata.tsv must contain pseudobulk_id")
    counts = counts.set_index("pseudobulk_id")
    meta = meta.set_index("pseudobulk_id")
    common = counts.index.intersection(meta.index)
    counts = counts.loc[common]
    meta = meta.loc[common]
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0)
    return counts, meta


def run_pseudobulk_logcpm_welch(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    condition_col: str,
    ctrl_group: str,
    test_group: str,
    min_count: int = 10,
    min_total_count: int = 15,
    min_samples_per_group: int = 2,
    prior_count: float = 0.5,
) -> PseudobulkDEResult:
    """Run a Python-native pseudobulk DE test on logCPM values.

    This intentionally keeps SCOOP's sample-level pseudobulk design. It is not a
    per-cell test and does not use integrated expression.
    """

    start = perf_counter()
    if condition_col not in metadata.columns:
        raise KeyError(f"condition_col {condition_col!r} is not present in metadata")
    keep_meta = metadata[metadata[condition_col].astype(str).isin([str(ctrl_group), str(test_group)])].copy()
    counts = counts.loc[keep_meta.index]
    condition = keep_meta[condition_col].astype(str)
    ctrl_mask = condition == str(ctrl_group)
    test_mask = condition == str(test_group)
    if int(ctrl_mask.sum()) < min_samples_per_group or int(test_mask.sum()) < min_samples_per_group:
        raise ValueError("not enough pseudobulk samples in one or both groups")

    values = counts.to_numpy(dtype=np.float64, copy=False)
    keep = filter_by_expr_like(values, ctrl_mask.to_numpy(), test_mask.to_numpy(), min_count=min_count, min_total_count=min_total_count)
    values = values[:, keep]
    genes = counts.columns[keep].astype(str)
    lib_size = values.sum(axis=1)
    lib_size = np.where(lib_size > 0, lib_size, 1.0)
    cpm = (values + prior_count) / (lib_size[:, None] + prior_count * values.shape[1]) * 1_000_000.0
    logcpm = np.log2(cpm)

    ctrl = logcpm[ctrl_mask.to_numpy(), :]
    test = logcpm[test_mask.to_numpy(), :]
    ctrl_mean = ctrl.mean(axis=0)
    test_mean = test.mean(axis=0)
    logfc = test_mean - ctrl_mean
    logcpm_mean = logcpm.mean(axis=0)
    ctrl_var = ctrl.var(axis=0, ddof=1)
    test_var = test.var(axis=0, ddof=1)
    n_ctrl = ctrl.shape[0]
    n_test = test.shape[0]
    se2 = ctrl_var / n_ctrl + test_var / n_test
    se = np.sqrt(np.maximum(se2, np.finfo(float).tiny))
    t_stat = logfc / se
    df_num = se2**2
    df_den = (ctrl_var / n_ctrl) ** 2 / max(n_ctrl - 1, 1) + (test_var / n_test) ** 2 / max(n_test - 1, 1)
    df = np.divide(df_num, df_den, out=np.full_like(df_num, max(n_ctrl + n_test - 2, 1), dtype=float), where=df_den > 0)
    pvalue = 2.0 * stats.t.sf(np.abs(t_stat), df=np.maximum(df, 1.0))
    pvalue = np.where(np.isfinite(pvalue), pvalue, 1.0)
    fdr = benjamini_hochberg(pvalue)

    out = pd.DataFrame(
        {
            "gene": genes,
            "logFC": logfc,
            "logCPM": logcpm_mean,
            "t": t_stat,
            "df": df,
            "PValue": pvalue,
            "FDR": fdr,
        }
    ).sort_values(["PValue", "gene"], kind="mergesort").reset_index(drop=True)
    design = pd.DataFrame(
        {
            "intercept": 1.0,
            f"{condition_col}_{test_group}_vs_{ctrl_group}": test_mask.astype(float).to_numpy(),
        },
        index=keep_meta.index,
    )
    manifest = {
        "schema_version": "scoop.pseudobulk_de.python_logcpm_welch.v1",
        "method": "python_logcpm_welch",
        "condition_col": condition_col,
        "ctrl_group": str(ctrl_group),
        "test_group": str(test_group),
        "n_samples": int(len(keep_meta)),
        "n_ctrl": int(ctrl_mask.sum()),
        "n_test": int(test_mask.sum()),
        "n_input_genes": int(counts.shape[1]),
        "n_tested_genes": int(out.shape[0]),
        "prior_count": float(prior_count),
        "min_count": int(min_count),
        "min_total_count": int(min_total_count),
        "seconds": round(perf_counter() - start, 6),
    }
    return PseudobulkDEResult(table=out, design=design, manifest=manifest)


def filter_by_expr_like(
    values: np.ndarray,
    ctrl_mask: np.ndarray,
    test_mask: np.ndarray,
    *,
    min_count: int,
    min_total_count: int,
) -> np.ndarray:
    n_min = int(min(ctrl_mask.sum(), test_mask.sum()))
    n_min = max(n_min, 1)
    expressed = ((values[ctrl_mask, :] >= min_count).sum(axis=0) >= n_min) | ((values[test_mask, :] >= min_count).sum(axis=0) >= n_min)
    total = values.sum(axis=0) >= min_total_count
    return expressed & total


def benjamini_hochberg(pvalue: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalue, dtype=float)
    n = p.size
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out = np.empty_like(adjusted)
    out[order] = adjusted
    return out


def compare_de_tables(candidate: pd.DataFrame, reference: pd.DataFrame, *, top_n: int = 100) -> dict:
    cand = candidate.set_index("gene")
    ref = reference.set_index("gene")
    common = cand.index.intersection(ref.index)
    if len(common) == 0:
        raise ValueError("candidate and reference have no common genes")
    logfc_corr = float(cand.loc[common, "logFC"].corr(ref.loc[common, "logFC"], method="spearman"))
    p_col = "PValue" if "PValue" in ref.columns else "PValue"
    score_col = "PValue"
    cand_top = set(cand.sort_values(score_col).head(top_n).index)
    ref_top = set(ref.sort_values(p_col).head(top_n).index)
    overlap = len(cand_top & ref_top) / max(1, min(top_n, len(cand_top), len(ref_top)))
    sign_agree = (np.sign(cand.loc[common, "logFC"]) == np.sign(ref.loc[common, "logFC"])).mean()
    return {
        "n_common_genes": int(len(common)),
        "spearman_logfc": logfc_corr,
        "top_n": int(top_n),
        "top_n_overlap_fraction": float(overlap),
        "logfc_sign_agreement": float(sign_agree),
    }
