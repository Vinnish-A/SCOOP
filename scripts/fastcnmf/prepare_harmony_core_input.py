#!/usr/bin/env python
"""Prepare identical core arrays for Harmony 0.2 vs 2.0 adapter tests."""
from __future__ import annotations

import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from cnmf.preprocess import stdscale_quantile_celing
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
RAW_ROOT = Path(os.environ.get("SCOOP_SPATIAL_RAW_ROOT", ROOT / ".scoop_local/data/raw/spatial/gbm_lowres_visium"))
OUT_DIR = ROOT / "tmp/fastcnmf_harmony2"
CORE_NPZ = OUT_DIR / "harmony_core_input.npz"
MANIFEST = OUT_DIR / "harmony_core_input_manifest.json"


def read_sample(sample_dir: Path) -> ad.AnnData:
    sample_id = sample_dir.name
    adata = sc.read_10x_h5(sample_dir / "filtered_feature_bc_matrix.h5")
    adata.var_names_make_unique()
    adata.obs_names = [f"{sample_id}:{barcode}" for barcode in adata.obs_names]
    adata.obs["sample_id"] = sample_id
    return adata


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    samples = sorted(p for p in RAW_ROOT.glob("GBM_*") if p.is_dir())
    adatas = [read_sample(sample) for sample in samples]
    merged = ad.concat(adatas, join="inner", label=None, merge="same")
    merged.X = sparse.csr_matrix(merged.X)
    sc.pp.filter_genes(merged, min_cells=20)

    sc.pp.highly_variable_genes(merged, flavor="seurat_v3", n_top_genes=3000)
    anorm = sc.pp.normalize_total(merged, target_sum=1e4, copy=True)
    anorm = anorm[:, merged.var["highly_variable"]]
    stdscale_quantile_celing(anorm, max_value=None, quantile_thresh=0.9999)
    sc.pp.pca(anorm, use_highly_variable=True, zero_center=True)

    hvg = merged[:, merged.var["highly_variable"]].copy()
    stdscale_quantile_celing(hvg, max_value=None, quantile_thresh=0.9999)
    x = np.asarray(hvg.X.todense() if sparse.issparse(hvg.X) else hvg.X, dtype=np.float64)
    pca = np.asarray(anorm.obsm["X_pca"], dtype=np.float64)
    sample_id = hvg.obs["sample_id"].astype(str).to_numpy()
    obs_names = hvg.obs_names.astype(str).to_numpy()
    var_names = hvg.var_names.astype(str).to_numpy()

    np.savez_compressed(
        CORE_NPZ,
        X=x,
        pca=pca,
        sample_id=sample_id,
        obs_names=obs_names,
        var_names=var_names,
    )
    manifest = {
        "core_npz": str(CORE_NPZ),
        "samples": [sample.name for sample in samples],
        "n_obs": int(x.shape[0]),
        "n_vars": int(x.shape[1]),
        "n_pcs": int(pca.shape[1]),
        "harmony_vars": ["sample_id"],
        "lamb": 1,
        "theta": 1,
        "max_iter_harmony": 20,
        "quantile_thresh": 0.9999,
        "n_top_genes": 3000,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
