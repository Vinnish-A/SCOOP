from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .omicverse_common import require_omicverse, unsupported_backend_result


def run_omicverse_mixed_core(adata, cfg: Mapping[str, Any], run_root: str | Path) -> dict[str, Any]:
    ov = require_omicverse()
    if hasattr(ov, "settings") and hasattr(ov.settings, "cpu_gpu_mixed_init"):
        ov.settings.cpu_gpu_mixed_init()
    return unsupported_backend_result("omicverse_cpu_gpu_mixed", cfg, run_root)
