#!/usr/bin/env python
"""Integrate sample-batch-adjusted cNMF consensus outputs into H5AD."""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TMP = ROOT / "tmp/cnmf_spatial_batch_adjusted_benchmark"
INPUT = TMP / "input/gbm_lowres_visium_3samples_sample_batch_adjusted.h5ad"
CNMF_OUT = TMP / "parallel/gbm_lowres_sample_batch_adjusted_cnmf"
OUTPUT = TMP / "integrated/gbm_lowres_visium_3samples_sample_batch_adjusted_with_cnmf.h5ad"
MANIFEST = TMP / "integrated/gbm_lowres_visium_3samples_sample_batch_adjusted_with_cnmf_manifest.json"


def main() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(INPUT)

    integrated = {
        "source_h5ad": str(INPUT),
        "cnmf_output": str(CNMF_OUT),
        "batch_adjustment": adata.uns.get("sample_batch_adjustment", {}),
        "k_values": {},
    }
    for k in (6, 8):
        usage_path = CNMF_OUT / f"gbm_lowres_sample_batch_adjusted_cnmf.usages.k_{k}.dt_0_5.consensus.txt"
        spectra_path = CNMF_OUT / f"gbm_lowres_sample_batch_adjusted_cnmf.spectra.k_{k}.dt_0_5.consensus.txt"

        usage = pd.read_csv(usage_path, sep="\t", index_col=0).reindex(adata.obs_names)
        usage.columns = [f"cnmf_batch_adjusted_k{k}_usage_{col}" for col in usage.columns]
        adata.obsm[f"X_cnmf_batch_adjusted_usage_k{k}"] = usage.to_numpy()

        spectra = pd.read_csv(spectra_path, sep="\t", index_col=0).reindex(columns=adata.var_names)
        spectra.index = [f"cnmf_batch_adjusted_k{k}_program_{idx}" for idx in spectra.index]
        adata.varm[f"cnmf_batch_adjusted_spectra_k{k}"] = spectra.T.to_numpy()
        adata.uns[f"cnmf_batch_adjusted_k{k}_usage_columns"] = list(usage.columns)
        adata.uns[f"cnmf_batch_adjusted_k{k}_spectra_columns"] = list(spectra.index)

        integrated["k_values"][str(k)] = {
            "usage": str(usage_path),
            "spectra": str(spectra_path),
            "obsm": f"X_cnmf_batch_adjusted_usage_k{k}",
            "varm": f"cnmf_batch_adjusted_spectra_k{k}",
            "usage_shape": list(usage.shape),
            "spectra_shape": list(spectra.shape),
        }

    adata.uns["cnmf_sample_batch_adjusted_benchmark"] = integrated
    adata.write_h5ad(OUTPUT)
    def json_default(value):
        if isinstance(value, np.integer):
            return int(value)
        if isinstance(value, np.floating):
            return float(value)
        if isinstance(value, np.ndarray):
            return value.tolist()
        raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")

    MANIFEST.write_text(json.dumps(integrated, indent=2, default=json_default), encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT), "shape": adata.shape, **integrated}, indent=2, default=json_default))


if __name__ == "__main__":
    main()
