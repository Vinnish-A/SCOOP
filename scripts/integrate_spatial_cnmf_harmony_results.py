#!/usr/bin/env python
"""Integrate Harmony-corrected cNMF consensus outputs into H5AD."""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TMP = ROOT / "tmp/cnmf_spatial_harmony_benchmark"
INPUT_BASE = TMP / "input/gbm_lowres_visium_3samples_harmony"
INPUT = Path(str(INPUT_BASE) + ".Corrected.HVG.Varnorm.h5ad")
TP10K = Path(str(INPUT_BASE) + ".TP10K.h5ad")
HVG = Path(str(INPUT_BASE) + ".Corrected.HVGs.txt")
CNMF_OUT = TMP / "parallel/gbm_lowres_harmony_cnmf"
OUTPUT = TMP / "integrated/gbm_lowres_visium_3samples_harmony_with_cnmf.h5ad"
MANIFEST = TMP / "integrated/gbm_lowres_visium_3samples_harmony_with_cnmf_manifest.json"


def json_default(value):
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def main() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(INPUT)
    adata.uns["cnmf_harmony_inputs"] = {
        "corrected_hvg_varnorm_h5ad": str(INPUT),
        "tp10k_h5ad": str(TP10K),
        "hvg_txt": str(HVG),
        "cnmf_output": str(CNMF_OUT),
    }

    integrated = {"source_h5ad": str(INPUT), "cnmf_output": str(CNMF_OUT), "k_values": {}}
    for k in (6, 8):
        usage_path = CNMF_OUT / f"gbm_lowres_harmony_cnmf.usages.k_{k}.dt_0_5.consensus.txt"
        spectra_path = CNMF_OUT / f"gbm_lowres_harmony_cnmf.spectra.k_{k}.dt_0_5.consensus.txt"

        usage = pd.read_csv(usage_path, sep="\t", index_col=0).reindex(adata.obs_names)
        usage.columns = [f"cnmf_harmony_k{k}_usage_{col}" for col in usage.columns]
        adata.obsm[f"X_cnmf_harmony_usage_k{k}"] = usage.to_numpy()

        spectra = pd.read_csv(spectra_path, sep="\t", index_col=0).reindex(columns=adata.var_names)
        spectra.index = [f"cnmf_harmony_k{k}_program_{idx}" for idx in spectra.index]
        adata.varm[f"cnmf_harmony_spectra_k{k}"] = spectra.T.to_numpy()
        adata.uns[f"cnmf_harmony_k{k}_usage_columns"] = list(usage.columns)
        adata.uns[f"cnmf_harmony_k{k}_spectra_columns"] = list(spectra.index)

        integrated["k_values"][str(k)] = {
            "usage": str(usage_path),
            "spectra": str(spectra_path),
            "obsm": f"X_cnmf_harmony_usage_k{k}",
            "varm": f"cnmf_harmony_spectra_k{k}",
            "usage_shape": list(usage.shape),
            "spectra_shape": list(spectra.shape),
        }

    adata.uns["cnmf_harmony_benchmark"] = integrated
    adata.write_h5ad(OUTPUT)
    MANIFEST.write_text(json.dumps(integrated, indent=2, default=json_default), encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT), "shape": adata.shape, **integrated}, indent=2, default=json_default))


if __name__ == "__main__":
    main()
