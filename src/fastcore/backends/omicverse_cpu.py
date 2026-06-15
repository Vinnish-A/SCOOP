from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from fastcore.vendor.omicverse_gpl import run_vendored_omicverse_cpu_core


def run_omicverse_cpu_core(adata, cfg: Mapping[str, Any], run_root: str | Path) -> dict[str, Any]:
    return run_vendored_omicverse_cpu_core(adata, cfg, run_root)
