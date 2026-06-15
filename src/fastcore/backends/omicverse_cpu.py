from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .omicverse_common import require_omicverse, unsupported_backend_result


def run_omicverse_cpu_core(adata, cfg: Mapping[str, Any], run_root: str | Path) -> dict[str, Any]:
    require_omicverse()
    return unsupported_backend_result("omicverse_cpu", cfg, run_root)
