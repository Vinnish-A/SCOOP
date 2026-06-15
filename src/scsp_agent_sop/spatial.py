from __future__ import annotations

from pathlib import Path
import shlex
import subprocess
from typing import Mapping, Any

import pandas as pd


def choose_rctd_mode(spatial_unit_type: str, expected_cells_per_unit: float | None = None, default_low_res: str = "full") -> str:
    unit = (spatial_unit_type or "").lower()
    if unit in {"roi", "spot", "large_bin"}:
        return default_low_res
    if unit in {"bin", "visium_hd_bin", "stereoseq_bin"}:
        if expected_cells_per_unit is not None and expected_cells_per_unit <= 3:
            return "multi"
        return default_low_res
    if unit in {"cell", "segmented_cell"}:
        return "doublet"
    return default_low_res


def build_rctd_command(command_template: str, **kwargs) -> list[str]:
    cmd = command_template.format(**kwargs)
    return shlex.split(cmd)


def run_rctd_command(command: list[str], cwd: str | Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(command, cwd=cwd, check=True, text=True, capture_output=True)


def summarize_rctd_weights(weights: pd.DataFrame) -> pd.DataFrame:
    """Summarise spatial unit x cell type weight table.

    The input must contain one row per spatial unit and one column per cell type.
    """
    numeric = weights.select_dtypes("number")
    if numeric.empty:
        raise ValueError("RCTD weights table has no numeric cell-type columns")
    top_type = numeric.idxmax(axis=1)
    top_weight = numeric.max(axis=1)
    entropy = -(numeric.div(numeric.sum(axis=1).replace(0, 1), axis=0) * (numeric.div(numeric.sum(axis=1).replace(0, 1), axis=0) + 1e-12).applymap(lambda x: __import__('math').log(x))).sum(axis=1)
    out = pd.DataFrame({
        "spatial_unit_id": weights.index.astype(str),
        "rctd_dominant_type": top_type.to_numpy(),
        "rctd_top1_weight": top_weight.to_numpy(float),
        "rctd_entropy": entropy.to_numpy(float),
        "rctd_n_types_weight_gt_0_05": (numeric > 0.05).sum(axis=1).to_numpy(int),
    })
    return out
