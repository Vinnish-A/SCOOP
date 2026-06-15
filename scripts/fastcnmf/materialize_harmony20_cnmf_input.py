#!/usr/bin/env python
"""Materialize harmonypy 2.0 adapter output as cNMF-compatible input files."""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TMP = ROOT / "tmp/fastcnmf_harmony2"
CORE_NPZ = TMP / "harmony_core_input.npz"
H20_NPZ = TMP / "harmony20_adapter_output.npz"
OUT_DIR = TMP / "input"
BASE = OUT_DIR / "gbm_lowres_harmony20_fastcnmf"
CORRECTED = Path(str(BASE) + ".Corrected.HVG.Varnorm.h5ad")
TP10K = Path(str(BASE) + ".TP10K.h5ad")
HVG = Path(str(BASE) + ".Corrected.HVGs.txt")
MANIFEST = OUT_DIR / "gbm_lowres_harmony20_fastcnmf_manifest.json"

BASELINE_INPUT = ROOT / "tmp/cnmf_spatial_harmony_benchmark/input"
BASELINE_TP10K = BASELINE_INPUT / "gbm_lowres_visium_3samples_harmony.TP10K.h5ad"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    core = np.load(CORE_NPZ, allow_pickle=True)
    h20 = np.load(H20_NPZ, allow_pickle=False)

    obs_names = core["obs_names"].astype(str)
    var_names = core["var_names"].astype(str)
    sample_id = core["sample_id"].astype(str)

    corrected = ad.AnnData(
        X=h20["X_corr"],
        obs=pd.DataFrame({"sample_id": sample_id}, index=obs_names),
        var=pd.DataFrame(index=var_names),
        obsm={"X_pca_harmony": h20["X_pca_harmony"]},
    )
    corrected.uns["fastcnmf_harmony2"] = {
        "method": "harmonypy2_adapter_fixed_lambda_moe",
        "lamb": 1,
        "theta": 1,
        "source_core_npz": str(CORE_NPZ),
    }
    corrected.write_h5ad(CORRECTED)

    tp10k = ad.read_h5ad(BASELINE_TP10K)
    tp10k = tp10k[obs_names, :].copy()
    tp10k.write_h5ad(TP10K)
    HVG.write_text("\n".join(var_names) + "\n", encoding="utf-8")

    manifest = {
        "corrected_h5ad": str(CORRECTED),
        "tp10k_h5ad": str(TP10K),
        "hvg_txt": str(HVG),
        "n_obs": int(corrected.n_obs),
        "n_vars": int(corrected.n_vars),
        "samples": sorted(set(sample_id)),
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
