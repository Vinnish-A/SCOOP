from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from scsp_agent_sop.config import deep_get
from scsp_agent_sop.omicverse_facilities import require_omicverse

from .omicverse_common import (
    build_options,
    configure_omicverse,
    copy_counts_to_x,
    run_harmony_bridge,
    run_omicverse_graph_umap_leiden,
    run_omicverse_preprocess_scale_pca,
    summarize_result,
    write_core_tables,
)


def run_omicverse_mixed_core(adata, cfg: Mapping[str, Any], run_root: str | Path) -> dict[str, Any]:
    ov = require_omicverse()
    options = build_options(cfg, adata, backend="fastcore_mixed", mode="cpu-gpu-mixed")
    configure_omicverse(ov, options)
    timings: dict[str, float] = {}
    copy_counts_to_x(adata, str(deep_get(cfg, "qc.counts_layer", "counts")))
    run_omicverse_preprocess_scale_pca(ov, adata, options, timings)
    harmony_used = run_harmony_bridge(adata, options, timings)
    leiden_results = run_omicverse_graph_umap_leiden(ov, adata, options, timings, cfg=cfg, run_root=run_root)
    artifacts = write_core_tables(
        adata, run_root, method="fastcore_mixed", options=options, leiden_results=leiden_results
    )
    return summarize_result(
        adata,
        options=options,
        timings=timings,
        harmony_used=harmony_used,
        artifacts=artifacts,
        quality_basis="fastcore_mixed_adapter_smoke",
    )
