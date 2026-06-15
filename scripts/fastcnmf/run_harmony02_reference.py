#!/usr/bin/env python
"""Run harmonypy 0.2 reference MOE correction on prepared arrays."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from fastcnmf.harmony2_compat import moe_correct_ridge_fast


TMP = ROOT / "tmp/fastcnmf_harmony2"
CORE_NPZ = TMP / "harmony_core_input.npz"
OUT_NPZ = TMP / "harmony02_reference_output.npz"
OUT_JSON = TMP / "harmony02_reference_metrics.json"


def main() -> None:
    import harmonypy

    data = np.load(CORE_NPZ, allow_pickle=True)
    x = data["X"]
    pca = data["pca"]
    obs = pd.DataFrame({"sample_id": data["sample_id"].astype(str)})

    start = time.perf_counter()
    harmony_res = harmonypy.run_harmony(
        pca,
        obs,
        "sample_id",
        lamb=1,
        theta=1,
        max_iter_harmony=20,
        verbose=False,
        random_state=0,
    )
    harmony_seconds = time.perf_counter() - start

    phi_moe = np.asarray(harmony_res.Phi_moe, dtype=np.float64)
    if phi_moe.shape[0] == x.shape[0]:
        phi_moe = phi_moe.T
    r = np.asarray(harmony_res.R, dtype=np.float64)
    if r.shape[0] != x.shape[0]:
        r = r.T
    start = time.perf_counter()
    x_corr = moe_correct_ridge_fast(x, r, phi_moe, np.asarray(harmony_res.lamb, dtype=np.float64))
    moe_seconds = time.perf_counter() - start
    z_corr = np.asarray(harmony_res.Z_corr, dtype=np.float64)
    if z_corr.shape[0] != x.shape[0]:
        z_corr = z_corr.T

    np.savez_compressed(OUT_NPZ, X_corr=x_corr, X_pca_harmony=z_corr, R=r)
    metrics = {
        "harmonypy_version": getattr(harmonypy, "__version__", "unknown"),
        "harmony_seconds": harmony_seconds,
        "moe_seconds": moe_seconds,
        "total_seconds": harmony_seconds + moe_seconds,
        "x_corr_shape": list(x_corr.shape),
        "r_shape": list(r.shape),
        "k": int(harmony_res.K),
    }
    OUT_JSON.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
