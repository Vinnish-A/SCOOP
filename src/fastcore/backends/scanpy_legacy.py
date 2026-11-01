from __future__ import annotations

from pathlib import Path
from time import perf_counter
from typing import Any, Mapping

from scsp_agent_sop.config import deep_get, resolve_obs_keys
from scsp_agent_sop.core import (
    build_identity_hvg_from_program_decision,
    core_batch_keys,
    correct_identity_pca,
    leiden_sweep,
    leiden_search_options,
    neighbors_umap,
    normalize_log1p,
    robust_z_by_sample,
    run_pca,
    score_programs,
    select_hvg,
)
from scsp_agent_sop.measurement import FeatureMeasurement, feature_measurement_path
from scsp_agent_sop.storage import register_file, write_table


def _time_step(timings: dict[str, float], name: str, fn, *args, **kwargs):
    start = perf_counter()
    out = fn(*args, **kwargs)
    timings[name] = perf_counter() - start
    return out


def run_scanpy_legacy_core(adata, cfg: Mapping[str, Any], run_root: str | Path) -> dict[str, Any]:
    """Run the original SCOOP core implementation as the single legacy fallback."""
    run_root = Path(run_root)
    timings: dict[str, float] = {}
    counts_layer = deep_get(cfg, "qc.counts_layer", "counts")
    sample_key = resolve_obs_keys(cfg, adata.obs).sample
    measurement_path = feature_measurement_path(adata)
    measurement = FeatureMeasurement.load(measurement_path) if measurement_path else None
    if measurement is not None:
        measurement.check(adata)

    _time_step(
        timings,
        "normalize_log1p",
        normalize_log1p,
        adata,
        counts_layer=counts_layer,
        target_sum=deep_get(cfg, "core.normalize_total_target_sum", 10000),
        size_factor_key=deep_get(cfg, "core.size_factor_key", None),
    )
    _time_step(
        timings,
        "score_programs",
        score_programs,
        adata,
        organism=deep_get(cfg, "run.organism", "human"),
        layer="log1p_norm",
        measurement=measurement,
    )
    _time_step(
        timings,
        "robust_z_by_sample",
        robust_z_by_sample,
        adata,
        ["stress_score", "ribo_score", "proliferation_score"],
        sample_key=sample_key,
    )

    hvg = _time_step(
        timings,
        "select_hvg_biology",
        select_hvg,
        adata,
        counts_layer=counts_layer,
        batch_key=sample_key,
        flavor=deep_get(cfg, "core.hvg_flavor", "seurat_v3"),
        n_top_genes=deep_get(cfg, "core.n_top_hvg", 3000),
        output_key="highly_variable_biology",
        span=float(deep_get(cfg, "core.hvg_span", 0.3)),
        retry_spans=deep_get(cfg, "core.hvg_retry_spans", []),
        fallback=deep_get(cfg, "core.hvg_fallback", "seurat"),
        min_cells=int(deep_get(cfg, "core.hvg_min_cells", 1)),
        measurement=measurement,
    )
    hvg_path = write_table(hvg, run_root / "core" / "tables" / "hvg_biology.parquet")
    register_file(adata, key="hvg_biology", path=hvg_path, schema="hvg_table.v1")

    _time_step(
        timings,
        "pca_biology",
        run_pca,
        adata,
        hvg_key="highly_variable_biology",
        obsm_key="X_pca_biology",
        n_comps=deep_get(cfg, "core.n_pcs", 50),
    )
    identity_hvg = _time_step(timings, "build_identity_hvg", build_identity_hvg_from_program_decision, adata)
    id_path = write_table(identity_hvg, run_root / "core" / "tables" / "hvg_identity.parquet")
    register_file(adata, key="hvg_identity", path=id_path, schema="hvg_identity.v1")
    _time_step(
        timings,
        "pca_identity_prebatch",
        run_pca,
        adata,
        hvg_key="highly_variable_identity",
        obsm_key="X_pca_identity_prebatch",
        n_comps=deep_get(cfg, "core.n_pcs", 50),
    )

    batch_keys, key_source = core_batch_keys(cfg, adata.obs)
    batch_method = str(deep_get(cfg, "core.batch_correction.method", "harmony2"))
    batch_correction = _time_step(
        timings,
        "harmony2",
        correct_identity_pca,
        adata,
        method=batch_method,
        batch_keys=batch_keys,
        key_source=key_source,
        max_iter_harmony=deep_get(cfg, "core.batch_correction.max_iter_harmony", 20),
        random_state=deep_get(cfg, "run.random_seed", 0),
        ncores=deep_get(cfg, "core.batch_correction.ncores", 0),
    )
    neighbors_n_pcs = deep_get(cfg, "core.neighbors_n_pcs", None)

    _time_step(
        timings,
        "neighbors_umap_identity",
        neighbors_umap,
        adata,
        use_rep="X_pca_harmony_identity",
        prefix="identity",
        n_neighbors=deep_get(cfg, "core.neighbors_n_neighbors", 15),
        n_pcs=neighbors_n_pcs,
        min_dist=deep_get(cfg, "core.umap_min_dist", 0.3),
        random_state=deep_get(cfg, "run.random_seed", 0),
    )
    _time_step(
        timings,
        "neighbors_umap_biology",
        neighbors_umap,
        adata,
        use_rep="X_pca_biology",
        prefix="biology",
        n_neighbors=deep_get(cfg, "core.neighbors_n_neighbors", 15),
        n_pcs=neighbors_n_pcs,
        min_dist=deep_get(cfg, "core.umap_min_dist", 0.3),
        random_state=deep_get(cfg, "run.random_seed", 0),
    )
    sweep, stability = _time_step(
        timings,
        "leiden_sweep",
        leiden_sweep,
        adata,
        graph_prefix="identity",
        resolutions=deep_get(cfg, "core.leiden_resolutions", [0.25, 0.5, 0.75, 1.0, 1.25, 1.5]),
        seeds=deep_get(cfg, "core.leiden_seeds", [0, 1, 2, 3, 4]),
        search_config=leiden_search_options(cfg, run_root),
    )
    sweep_path = write_table(sweep, run_root / "core" / "tables" / "leiden_sweep.parquet")
    st_path = write_table(stability, run_root / "core" / "tables" / "cluster_stability.parquet")
    register_file(adata, key="leiden_sweep", path=sweep_path, schema="leiden_sweep.v1")
    register_file(adata, key="cluster_stability", path=st_path, schema="cluster_stability.v1")

    return {
        "backend": "scanpy_legacy",
        "harmony2_used": batch_correction["status"] == "corrected",
        "batch_correction_method": batch_method,
        "batch_keys": batch_correction["keys"],
        "batch_correction": batch_correction,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_clusters": int(adata.obs["cluster_identity"].nunique()) if "cluster_identity" in adata.obs else None,
        "timings": timings,
        "artifacts": {
            "hvg_biology": str(hvg_path),
            "hvg_identity": str(id_path),
            "leiden_sweep": str(sweep_path),
            "cluster_stability": str(st_path),
            **({"feature_measurement": str(measurement_path)} if measurement_path else {}),
        },
    }
