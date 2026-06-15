from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .omicverse_common import (
    build_options,
    configure_omicverse,
    map_standard_core_keys,
    require_omicverse,
    run_harmony_bridge,
    run_omicverse_graph_umap_leiden,
    run_omicverse_preprocess_scale_pca,
    summarize_result,
    time_step,
    write_core_tables,
)
from scsp_agent_sop.config import deep_get
from scsp_agent_sop.storage import init_file_registry


def run_omicverse_rust_oom_core(input_h5ad: str | Path, output_h5ad: str | Path, cfg: Mapping[str, Any], run_root: str | Path) -> dict[str, Any]:
    ov = require_omicverse()
    if not hasattr(ov, "read"):
        raise RuntimeError("OmicVerse rust/OOM backend requires ov.read(path, backend='rust').")
    run_root = Path(run_root)
    output_h5ad = Path(output_h5ad)
    timings: dict[str, float] = {}
    adata_oom = time_step(timings, "ov_read_rust", ov.read, str(input_h5ad), backend="rust")
    options = build_options(cfg, adata_oom, backend="omicverse_rust_oom", mode="cpu")
    rust_mode = str(deep_get(cfg, "core.fastcore.omicverse.rust_oom.preprocess_mode", "shiftlog|pearson"))
    if "|" not in rust_mode:
        rust_mode = f"{rust_mode}|pearson"
    options = type(options)(**{**options.__dict__, "preprocess_mode": rust_mode, "use_implicit_centering": False})
    configure_omicverse(ov, options)
    try:
        run_omicverse_preprocess_scale_pca(ov, adata_oom, options, timings)
        materialize = getattr(adata_oom, "to_adata", None) or getattr(adata_oom, "to_memory", None)
        if materialize is None:
            raise RuntimeError("OmicVerse Rust/OOM object must expose to_adata() or to_memory().")
        adata = time_step(timings, "materialize_minimal_adata", materialize)
        init_file_registry(adata, str(deep_get(cfg, "run.run_id", run_root.name)))
        map_standard_core_keys(adata)
        harmony_used = run_harmony_bridge(adata, options, timings)
        run_omicverse_graph_umap_leiden(ov, adata, options, timings)
        artifacts = write_core_tables(adata, run_root, method="omicverse_rust_oom", options=options)
        result = summarize_result(
            adata,
            options=options,
            timings=timings,
            harmony_used=harmony_used,
            artifacts=artifacts,
            quality_basis="omicverse_rust_oom_adapter_smoke",
        )
        result["_adata"] = adata
        result["output_h5ad"] = str(output_h5ad)
        return result
    finally:
        close = getattr(adata_oom, "close", None)
        if close is not None:
            close()
