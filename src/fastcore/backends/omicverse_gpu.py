from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .omicverse_common import require_omicverse, unsupported_backend_result


def run_omicverse_gpu_core(adata, cfg: Mapping[str, Any], run_root: str | Path) -> dict[str, Any]:
    ov = require_omicverse()
    moved = False
    try:
        if hasattr(ov, "pp") and hasattr(ov.pp, "anndata_to_GPU"):
            ov.pp.anndata_to_GPU(adata)
            moved = True
        return unsupported_backend_result("omicverse_gpu_rapids", cfg, run_root)
    finally:
        if moved and hasattr(ov, "pp") and hasattr(ov.pp, "anndata_to_CPU"):
            ov.pp.anndata_to_CPU(adata)
