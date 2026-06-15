#!/usr/bin/env python
"""Build a combined low-resolution Visium H5AD for cNMF benchmarking."""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
RAW_ROOT = ROOT / "data/raw/spatial/gbm_lowres_visium"
OUT_DIR = ROOT / "tmp/cnmf_spatial_benchmark/input"
OUT_H5AD = OUT_DIR / "gbm_lowres_visium_3samples_cnmf_input.h5ad"
MANIFEST = OUT_DIR / "gbm_lowres_visium_3samples_cnmf_input_manifest.json"


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

    # Keep genes with enough support for a quick but realistic cNMF benchmark.
    sc.pp.filter_genes(combined, min_cells=20)
    sc.pp.highly_variable_genes(
        combined,
        flavor="seurat_v3",
        n_top_genes=min(3000, combined.n_vars),
        layer="counts",
        batch_key="sample_id",
    )
    combined = combined[:, combined.var["highly_variable"].to_numpy()].copy()
    combined.layers["counts"] = combined.X.copy()
    combined.uns["schema_version"] = "scsp_agent_sop.h5ad.v1"

    combined.write_h5ad(OUT_H5AD)
    summary = {
        "input_root": str(RAW_ROOT),
        "output_h5ad": str(OUT_H5AD),
        "samples": [sample.name for sample in samples],
        "n_obs": int(combined.n_obs),
        "n_vars": int(combined.n_vars),
        "spots_by_sample": combined.obs["sample_id"].value_counts().sort_index().to_dict(),
        "layers": list(combined.layers.keys()),
    }
    MANIFEST.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
