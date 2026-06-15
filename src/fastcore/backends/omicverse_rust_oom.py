from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .omicverse_common import require_omicverse, unsupported_backend_result


def run_omicverse_rust_oom_core(input_h5ad: str | Path, output_h5ad: str | Path, cfg: Mapping[str, Any], run_root: str | Path) -> dict[str, Any]:
    ov = require_omicverse()
    if not hasattr(ov, "read"):
        raise RuntimeError("OmicVerse rust/OOM backend requires ov.read(path, backend='rust').")
    return unsupported_backend_result("omicverse_rust_oom", cfg, run_root)
