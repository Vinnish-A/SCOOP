#!/usr/bin/env python
"""Run harmonypy 2.0 with FastCNMF's cNMF MOE compatibility adapter."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from fastcnmf.harmony2_compat import harmony2_moe_correct


TMP = ROOT / "tmp/fastcnmf_harmony2"
CORE_NPZ = TMP / "harmony_core_input.npz"
OUT_NPZ = TMP / "harmony20_adapter_output.npz"
OUT_JSON = TMP / "harmony20_adapter_metrics.json"


def main() -> None:
    import harmonypy

    data = np.load(CORE_NPZ, allow_pickle=True)
    x = data["X"]
    pca = data["pca"]
    obs = pd.DataFrame({"sample_id": data["sample_id"].astype(str)})

    start = time.perf_counter()
    result = harmony2_moe_correct(
        x,
        pca,
        obs,
        "sample_id",
        lamb=1,
        theta=1,
        max_iter_harmony=20,
        ncores=1,
        verbose=False,
        random_state=0,
    )
    total_seconds = time.perf_counter() - start
    np.savez_compressed(
        OUT_NPZ,
        X_corr=result.x_corr,
        X_pca_harmony=result.x_pca_harmony,
        R=result.r,
        Phi_moe=result.phi_moe,
        lamb=result.lamb,
    )
    metrics = {
        "harmonypy_version": getattr(harmonypy, "__version__", "unknown"),
        "total_seconds": total_seconds,
        "x_corr_shape": list(result.x_corr.shape),
        "r_shape": list(result.r.shape),
        "phi_moe_shape": list(result.phi_moe.shape),
        "k": int(result.k),
    }
    OUT_JSON.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
