#!/usr/bin/env python
"""Build a merged Visium H5AD with sample-level batch scaling removed.

The output X matrix is nonnegative and suitable for cNMF. For each gene, counts
within each sample are scaled so the sample mean matches the global mean.
"""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
RAW_ROOT = ROOT / "data/raw/spatial/gbm_lowres_visium"
OUT_DIR = ROOT / "tmp/cnmf_spatial_batch_adjusted_benchmark/input"
OUT_H5AD = OUT_DIR / "gbm_lowres_visium_3samples_sample_batch_adjusted.h5ad"
MANIFEST = OUT_DIR / "gbm_lowres_visium_3samples_sample_batch_adjusted_manifest.json"


def read_sample(sample_dir: Path) -> ad.AnnData:
    sample_id = sample_dir.name
    adata = sc.read_10x_h5(sample_dir / "filtered_feature_bc_matrix.h5")
    adata.var_names_make_unique()
    adata.obs_names = [f"{sample_id}:{barcode}" for barcode in adata.obs_names]
    adata.obs["sample_id"] = sample_id
    adata.obs["spatial_unit_type"] = "spot"

    positions = pd.read_csv(sample_dir / "spatial/tissue_positions.csv")
    barcode_col = positions.columns[0]
    positions.index = [f"{sample_id}:{barcode}" for barcode in positions[barcode_col].astype(str)]
    positions = positions.reindex(adata.obs_names)
    adata.obs["in_tissue"] = positions["in_tissue"].astype("Int64").astype(str).values
    adata.obs["array_row"] = positions["array_row"].values
    adata.obs["array_col"] = positions["array_col"].values
    adata.obs["pxl_row_in_fullres"] = positions["pxl_row_in_fullres"].values
    adata.obs["pxl_col_in_fullres"] = positions["pxl_col_in_fullres"].values
    adata.obsm["spatial"] = positions[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
    return adata


def remove_sample_batch_by_gene_scaling(adata: ad.AnnData) -> ad.AnnData:
    raw = sparse.csr_matrix(adata.X)
    global_mean = np.asarray(raw.mean(axis=0)).ravel()
    adjusted_parts = []
    sample_summaries = {}

    for sample_id in adata.obs["sample_id"].astype(str).unique():
        mask = adata.obs["sample_id"].astype(str).to_numpy() == sample_id
        sample_matrix = raw[mask].tocsr()
        sample_mean = np.asarray(sample_matrix.mean(axis=0)).ravel()
        ratio = np.ones_like(global_mean, dtype=np.float64)
        nonzero = sample_mean > 0
        ratio[nonzero] = global_mean[nonzero] / sample_mean[nonzero]
        adjusted = sample_matrix.multiply(ratio).tocsr()
        adjusted_parts.append(adjusted)
        sample_summaries[sample_id] = {
            "spots": int(mask.sum()),
            "genes_with_nonzero_sample_mean": int(nonzero.sum()),
            "median_nonzero_scaling_ratio": float(np.median(ratio[nonzero])),
        }

    out = adata.copy()
    out.layers["raw_counts"] = raw.copy()
    out.X = sparse.vstack(adjusted_parts, format="csr")
    out.layers["sample_batch_adjusted"] = out.X.copy()
    out.uns["sample_batch_adjustment"] = {
        "method": "per_gene_sample_mean_scaling_to_global_mean",
        "batch_key": "sample_id",
        "matrix": "X",
        "raw_counts_layer": "raw_counts",
        "summaries": sample_summaries,
    }
    return out


def select_variable_genes(adata: ad.AnnData, n_top: int = 3000) -> ad.AnnData:
    norm = adata.copy()
    sc.pp.normalize_total(norm, target_sum=1e4)
    sc.pp.log1p(norm)
    x = norm.X.toarray() if sparse.issparse(norm.X) else np.asarray(norm.X)
    variances = np.var(x, axis=0)
    top = np.argsort(variances)[::-1][: min(n_top, adata.n_vars)]
    keep = np.zeros(adata.n_vars, dtype=bool)
    keep[top] = True
    adata.var["batch_adjusted_log1p_variance"] = variances
    adata.var["highly_variable"] = keep
    return adata[:, keep].copy()


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    samples = sorted(p for p in RAW_ROOT.glob("GBM_*") if p.is_dir())
    if not samples:
        raise FileNotFoundError(f"No GBM sample directories found under {RAW_ROOT}")

    adatas = [read_sample(sample) for sample in samples]
    combined = ad.concat(adatas, join="inner", label=None, merge="same")
    combined.X = sparse.csr_matrix(combined.X)
    combined.layers["counts"] = combined.X.copy()
    combined.obs["total_counts"] = combined.X.sum(axis=1).A1
    combined.obs["n_genes_by_counts"] = (combined.X > 0).sum(axis=1).A1

    sc.pp.filter_genes(combined, min_cells=20)
    adjusted = remove_sample_batch_by_gene_scaling(combined)
    adjusted = select_variable_genes(adjusted, n_top=3000)
    adjusted.layers["sample_batch_adjusted"] = adjusted.X.copy()
    adjusted.uns["schema_version"] = "scsp_agent_sop.h5ad.v1"

    adjusted.write_h5ad(OUT_H5AD)
    summary = {
        "input_root": str(RAW_ROOT),
        "output_h5ad": str(OUT_H5AD),
        "samples": [sample.name for sample in samples],
        "n_obs": int(adjusted.n_obs),
        "n_vars": int(adjusted.n_vars),
        "spots_by_sample": adjusted.obs["sample_id"].value_counts().sort_index().to_dict(),
        "batch_key_removed": "sample_id",
        "batch_adjustment": adjusted.uns["sample_batch_adjustment"],
        "x_matrix": "sample_batch_adjusted_nonnegative",
        "layers": list(adjusted.layers.keys()),
    }
    MANIFEST.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
