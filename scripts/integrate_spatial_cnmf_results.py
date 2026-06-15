#!/usr/bin/env python
"""Integrate cNMF consensus outputs into the benchmark spatial H5AD."""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TMP = ROOT / "tmp/cnmf_spatial_benchmark"
INPUT = TMP / "input/gbm_lowres_visium_3samples_cnmf_input.h5ad"
CNMF_OUT = TMP / "parallel/gbm_lowres_cnmf"
OUTPUT = TMP / "integrated/gbm_lowres_visium_3samples_with_parallel_cnmf.h5ad"
MANIFEST = TMP / "integrated/gbm_lowres_visium_3samples_with_parallel_cnmf_manifest.json"


def main() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(INPUT)

    integrated = {"source_h5ad": str(INPUT), "cnmf_output": str(CNMF_OUT), "k_values": {}}
    for k in (6, 8):
        usage_path = CNMF_OUT / f"gbm_lowres_cnmf.usages.k_{k}.dt_0_5.consensus.txt"
        spectra_path = CNMF_OUT / f"gbm_lowres_cnmf.spectra.k_{k}.dt_0_5.consensus.txt"

        usage = pd.read_csv(usage_path, sep="\t", index_col=0)
        usage = usage.reindex(adata.obs_names)
        usage.columns = [f"cnmf_k{k}_usage_{col}" for col in usage.columns]
        adata.obsm[f"X_cnmf_usage_k{k}"] = usage.to_numpy()

        spectra = pd.read_csv(spectra_path, sep="\t", index_col=0)
        spectra = spectra.reindex(columns=adata.var_names)
        spectra.index = [f"cnmf_k{k}_program_{idx}" for idx in spectra.index]
        adata.varm[f"cnmf_spectra_k{k}"] = spectra.T.to_numpy()
        adata.uns[f"cnmf_k{k}_usage_columns"] = list(usage.columns)
        adata.uns[f"cnmf_k{k}_spectra_columns"] = list(spectra.index)

        integrated["k_values"][str(k)] = {
            "usage": str(usage_path),
            "spectra": str(spectra_path),
            "obsm": f"X_cnmf_usage_k{k}",
            "varm": f"cnmf_spectra_k{k}",
            "usage_shape": list(usage.shape),
            "spectra_shape": list(spectra.shape),
        }

    adata.uns["cnmf_benchmark"] = integrated
    adata.write_h5ad(OUTPUT)
    MANIFEST.write_text(json.dumps(integrated, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT), "shape": adata.shape, **integrated}, indent=2))


if __name__ == "__main__":
    main()
