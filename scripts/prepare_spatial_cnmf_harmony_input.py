#!/usr/bin/env python
"""Prepare merged Visium input using cNMF's Harmony preprocessing."""
from __future__ import annotations

import json
import os
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
from cnmf.preprocess import Preprocess
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
RAW_ROOT = Path(os.environ.get("SCOOP_SPATIAL_RAW_ROOT", ROOT / ".scoop_local/data/raw/spatial/gbm_lowres_visium"))
OUT_DIR = ROOT / "tmp/cnmf_spatial_harmony_benchmark/input"
MERGED_H5AD = OUT_DIR / "gbm_lowres_visium_3samples_merged_raw.h5ad"
BASE = OUT_DIR / "gbm_lowres_visium_3samples_harmony"
CORRECTED_H5AD = Path(str(BASE) + ".Corrected.HVG.Varnorm.h5ad")
TP10K_H5AD = Path(str(BASE) + ".TP10K.h5ad")
HVG_TXT = Path(str(BASE) + ".Corrected.HVGs.txt")
MANIFEST = OUT_DIR / "gbm_lowres_visium_3samples_harmony_manifest.json"


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
    adatas = [read_sample(sample) for sample in samples]
    merged = ad.concat(adatas, join="inner", label=None, merge="same")
    merged.X = sparse.csr_matrix(merged.X)
    merged.layers["counts"] = merged.X.copy()
    merged.obs["total_counts"] = merged.X.sum(axis=1).A1
    merged.obs["n_genes_by_counts"] = (merged.X > 0).sum(axis=1).A1
    sc.pp.filter_genes(merged, min_cells=20)
    merged.uns["schema_version"] = "scsp_agent_sop.h5ad.v1"
    merged.write_h5ad(MERGED_H5AD)

    pre = Preprocess(random_seed=20260614)
    corrected, tp10k, hvgs = pre.preprocess_for_cnmf(
        merged,
        harmony_vars="sample_id",
        n_top_rna_genes=3000,
        librarysize_targetsum=1e4,
        makeplots=False,
        theta=1,
        save_output_base=str(BASE),
        max_iter_harmony=20,
    )
    corrected.uns["schema_version"] = "scsp_agent_sop.h5ad.v1"
    corrected.uns["cnmf_harmony_preprocess"] = {
        "method": "cnmf.preprocess.Preprocess.preprocess_for_cnmf",
        "harmony_vars": ["sample_id"],
        "n_top_rna_genes": 3000,
        "librarysize_targetsum": 1e4,
        "theta": 1,
        "max_iter_harmony": 20,
        "tp10k_h5ad": str(TP10K_H5AD),
        "hvg_txt": str(HVG_TXT),
    }
    corrected.write_h5ad(CORRECTED_H5AD)

    summary = {
        "input_root": str(RAW_ROOT),
        "merged_raw_h5ad": str(MERGED_H5AD),
        "corrected_h5ad": str(CORRECTED_H5AD),
        "tp10k_h5ad": str(TP10K_H5AD),
        "hvg_txt": str(HVG_TXT),
        "samples": [sample.name for sample in samples],
        "n_obs": int(corrected.n_obs),
        "n_vars": int(corrected.n_vars),
        "spots_by_sample": corrected.obs["sample_id"].value_counts().sort_index().to_dict(),
        "batch_key_removed": "sample_id",
        "method": "cNMF Preprocess Harmony MOE ridge correction",
        "obsm": list(corrected.obsm.keys()),
    }
    MANIFEST.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
