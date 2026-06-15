from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .omicverse_common import (
    build_options,
    configure_omicverse,
    copy_counts_to_x,
    require_omicverse,
    run_harmony_bridge,
    run_omicverse_graph_umap_leiden,
    run_omicverse_preprocess_scale_pca,
    summarize_result,
    time_step,
    write_core_tables,
)
from scsp_agent_sop.config import deep_get


def run_omicverse_gpu_core(adata, cfg: Mapping[str, Any], run_root: str | Path) -> dict[str, Any]:
    ov = require_omicverse()
    options = build_options(cfg, adata, backend="omicverse_gpu_rapids", mode="gpu")
    configure_omicverse(ov, options)
    timings: dict[str, float] = {}
    moved = False
    try:
        copy_counts_to_x(adata, str(deep_get(cfg, "qc.counts_layer", "counts")))
        time_step(timings, "anndata_to_gpu_preprocess", ov.pp.anndata_to_GPU, adata)
        moved = True
        run_omicverse_preprocess_scale_pca(ov, adata, options, timings)

        # Harmony 2.0 is CPU-only, so bridge through CPU after PCA, then move
        # the corrected embedding back to GPU for graph/UMAP/Leiden.
        time_step(timings, "anndata_to_cpu_harmony_bridge", ov.pp.anndata_to_CPU, adata)
        moved = False
        harmony_used = run_harmony_bridge(adata, options, timings)
        time_step(timings, "anndata_to_gpu_graph", ov.pp.anndata_to_GPU, adata)
        moved = True
        run_omicverse_graph_umap_leiden(ov, adata, options, timings)
        time_step(timings, "anndata_to_cpu_final", ov.pp.anndata_to_CPU, adata)
        moved = False
        artifacts = write_core_tables(adata, run_root, method="omicverse_gpu_rapids", options=options)
        return summarize_result(
            adata,
            options=options,
            timings=timings,
            harmony_used=harmony_used,
            artifacts=artifacts,
            quality_basis="omicverse_gpu_rapids_adapter_smoke",
        )
    finally:
        if moved and hasattr(ov, "pp") and hasattr(ov.pp, "anndata_to_CPU"):
            ov.pp.anndata_to_CPU(adata)
