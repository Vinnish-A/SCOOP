from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np
import pandas as pd

from .matrix import get_matrix, sum_axis, nnz_axis, subset_gene_sum


@dataclass
class ThresholdConfig:
    low_counts_mad: float = 3.0
    low_genes_mad: float = 3.0
    high_mt_mad: float = 3.0
    high_ribo_mad: float = 4.0
    high_hb_mad: float = 4.0
    q_low: float = 0.10
    q_high: float = 0.90
    q_high_strict: float = 0.95


def annotate_gene_flags(adata, organism: str = "human", symbol_key: str | None = None) -> None:
    symbols = adata.var[symbol_key].astype(str) if symbol_key and symbol_key in adata.var else adata.var_names.astype(str)
    if organism.lower().startswith("mouse"):
        mt = symbols.str.startswith("mt-")
        ribo = symbols.str.startswith(("Rpl", "Rps"))
        hb = symbols.str.match(r"^Hb[ab].*|^Hbb.*|^Hba.*", case=False)
    else:
        upper = symbols.str.upper()
        mt = upper.str.startswith("MT-")
        ribo = upper.str.startswith(("RPL", "RPS"))
        hb = upper.str.match(r"^HB[ABDEGQMZ].*")
    adata.var["mt_gene"] = np.asarray(mt, dtype=bool)
    adata.var["ribo_gene"] = np.asarray(ribo, dtype=bool)
    adata.var["hb_gene"] = np.asarray(hb, dtype=bool)


def compute_basic_qc(adata, counts_layer: str = "counts") -> None:
    X = get_matrix(adata, counts_layer)
    total = sum_axis(X, axis=1).astype(float)
    genes = nnz_axis(X, axis=1).astype(int)
    adata.obs["total_counts"] = total
    adata.obs["n_genes_by_counts"] = genes
    denom = np.maximum(total, 1.0)
    for flag, col in [
        ("mt_gene", "pct_counts_mt"),
        ("ribo_gene", "pct_counts_ribo"),
        ("hb_gene", "pct_counts_hb"),
    ]:
        if flag not in adata.var:
            adata.obs[col] = 0.0
        else:
            s = subset_gene_sum(adata, np.asarray(adata.var[flag], dtype=bool), counts_layer)
            adata.obs[col] = 100.0 * s / denom
    adata.obs["log10_total_counts"] = np.log10(total + 1)
    adata.obs["log10_n_genes_by_counts"] = np.log10(genes + 1)


def _mad(x: np.ndarray) -> float:
    med = np.nanmedian(x)
    return float(1.4826 * np.nanmedian(np.abs(x - med)))


def _thresholds(x: np.ndarray) -> dict[str, float]:
    med = float(np.nanmedian(x))
    mad = max(_mad(x), 1e-9)
    return {
        "median": med,
        "mad": mad,
        "q05": float(np.nanquantile(x, 0.05)),
        "q10": float(np.nanquantile(x, 0.10)),
        "q90": float(np.nanquantile(x, 0.90)),
        "q95": float(np.nanquantile(x, 0.95)),
    }


def assign_qc_flags_per_sample(
    adata,
    sample_key: str = "sample_id",
    cfg: ThresholdConfig | None = None,
) -> pd.DataFrame:
    cfg = cfg or ThresholdConfig()
    required = ["log10_total_counts", "log10_n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"]
    for col in required:
        if col not in adata.obs:
            raise KeyError(f"Missing QC column: {col}. Run compute_basic_qc first.")
    for col in [
        "qc_low_counts_flag", "qc_low_genes_flag", "qc_high_mt_flag",
        "qc_high_ribo_flag", "qc_high_hb_flag", "qc_pass", "qc_class",
        "qc_failure_reason",
    ]:
        adata.obs[col] = False if col.endswith("flag") or col == "qc_pass" else ""
    rows: list[dict[str, Any]] = []
    for sample, idx in adata.obs.groupby(sample_key).indices.items():
        idx_list = list(idx)
        obs = adata.obs.iloc[idx_list]
        stat = {col: _thresholds(obs[col].to_numpy(float)) for col in required}
        low_counts = (obs["log10_total_counts"] < stat["log10_total_counts"]["median"] - cfg.low_counts_mad * stat["log10_total_counts"]["mad"]) & (obs["log10_total_counts"] < stat["log10_total_counts"]["q10"])
        low_genes = (obs["log10_n_genes_by_counts"] < stat["log10_n_genes_by_counts"]["median"] - cfg.low_genes_mad * stat["log10_n_genes_by_counts"]["mad"]) & (obs["log10_n_genes_by_counts"] < stat["log10_n_genes_by_counts"]["q10"])
        high_mt = (obs["pct_counts_mt"] > stat["pct_counts_mt"]["median"] + cfg.high_mt_mad * stat["pct_counts_mt"]["mad"]) & (obs["pct_counts_mt"] > stat["pct_counts_mt"]["q90"])
        high_ribo = (obs["pct_counts_ribo"] > stat["pct_counts_ribo"]["median"] + cfg.high_ribo_mad * stat["pct_counts_ribo"]["mad"]) & (obs["pct_counts_ribo"] > stat["pct_counts_ribo"]["q95"])
        high_hb = (obs["pct_counts_hb"] > stat["pct_counts_hb"]["median"] + cfg.high_hb_mad * stat["pct_counts_hb"]["mad"]) & (obs["pct_counts_hb"] > stat["pct_counts_hb"]["q95"])
        adata.obs.loc[obs.index, "qc_low_counts_flag"] = low_counts.to_numpy(bool)
        adata.obs.loc[obs.index, "qc_low_genes_flag"] = low_genes.to_numpy(bool)
        adata.obs.loc[obs.index, "qc_high_mt_flag"] = high_mt.to_numpy(bool)
        adata.obs.loc[obs.index, "qc_high_ribo_flag"] = high_ribo.to_numpy(bool)
        adata.obs.loc[obs.index, "qc_high_hb_flag"] = high_hb.to_numpy(bool)
        fail = (low_counts & low_genes) | (high_mt & low_genes) | (high_mt & low_counts)
        suspect = (~fail) & (high_mt | high_ribo | high_hb | low_counts | low_genes)
        adata.obs.loc[obs.index, "qc_class"] = np.where(fail, "fail", np.where(suspect, "suspect", "pass"))
        adata.obs.loc[obs.index, "qc_pass"] = ~fail
        reasons = []
        for a, b, c, d, e in zip(low_counts, low_genes, high_mt, high_ribo, high_hb):
            r = []
            if a: r.append("low_counts")
            if b: r.append("low_genes")
            if c: r.append("high_mt")
            if d: r.append("high_ribo")
            if e: r.append("high_hb")
            reasons.append(";".join(r))
        adata.obs.loc[obs.index, "qc_failure_reason"] = reasons
        row = {"sample_id": sample, "n_obs": len(obs), "fail_fraction": float(fail.mean()), "suspect_fraction": float(suspect.mean())}
        for col, st in stat.items():
            for k, v in st.items():
                row[f"{col}_{k}"] = v
        rows.append(row)
    adata.obs["final_use"] = adata.obs["qc_pass"].astype(bool)
    return pd.DataFrame(rows)


def run_scrublet_per_sample(
    adata,
    sample_key: str = "sample_id",
    counts_layer: str = "counts",
    expected_doublet_rate: float = 0.05,
    stdev_doublet_rate: float = 0.02,
    sim_doublet_ratio: float = 2.0,
    random_state: int = 0,
) -> pd.DataFrame:
    import scanpy as sc

    adata.obs["doublet_score"] = np.nan
    adata.obs["doublet_call_scrublet"] = False
    adata.obs["doublet_call"] = "not_run"
    rows = []
    for sample, names in adata.obs.groupby(sample_key).groups.items():
        ad = adata[names].copy()
        ad.X = ad.layers[counts_layer].copy()
        sc.pp.scrublet(
            ad,
            expected_doublet_rate=expected_doublet_rate,
            stdev_doublet_rate=stdev_doublet_rate,
            sim_doublet_ratio=sim_doublet_ratio,
            random_state=random_state,
            threshold=None,
        )
        scores = ad.obs["doublet_score"].astype(float)
        pred = ad.obs["predicted_doublet"].astype(bool)
        adata.obs.loc[ad.obs_names, "doublet_score"] = scores.to_numpy()
        adata.obs.loc[ad.obs_names, "doublet_call_scrublet"] = pred.to_numpy()
        high_complexity = adata.obs.loc[ad.obs_names, "n_genes_by_counts"] > adata.obs.loc[ad.obs_names, "n_genes_by_counts"].quantile(0.95)
        high_counts = adata.obs.loc[ad.obs_names, "total_counts"] > adata.obs.loc[ad.obs_names, "total_counts"].quantile(0.95)
        confident = pred & (high_counts | high_complexity)
        call = np.where(confident, "doublet", np.where(pred, "ambiguous", "singlet"))
        adata.obs.loc[ad.obs_names, "doublet_call"] = call
        rows.append({
            "sample_id": sample,
            "n_obs": ad.n_obs,
            "scrublet_predicted_fraction": float(pred.mean()),
            "high_confidence_doublet_fraction": float(confident.mean()),
            "score_median": float(scores.median()),
            "score_q95": float(scores.quantile(0.95)),
        })
    # Remove only high-confidence doublets from default downstream.
    adata.obs.loc[adata.obs["doublet_call"] == "doublet", "final_use"] = False
    return pd.DataFrame(rows)
