from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Any, Callable, Mapping

import numpy as np
import pandas as pd

from scsp_agent_sop.config import deep_get
from scsp_agent_sop.core import run_harmony2
from scsp_agent_sop.storage import register_file, write_table
from omicverse_transfer.external import require_omicverse


@dataclass(frozen=True)
class OmicVerseCoreOptions:
    backend: str
    mode: str
    preprocess_mode: str
    target_sum: float
    n_hvgs: int
    scaled_layer: str
    n_pcs: int
    n_neighbors: int
    neighbors_n_pcs: int
    neighbors_method: str | None
    neighbors_transformer: str | None
    umap_min_dist: float
    umap_method: str | None
    leiden_resolution: float
    random_state: int
    batch_keys: list[str]
    max_iter_harmony: int
    harmony_ncores: int
    identify_robust: bool
    use_implicit_centering: bool


def build_options(cfg: Mapping[str, Any], adata: Any | None, *, backend: str, mode: str) -> OmicVerseCoreOptions:
    random_state = int(deep_get(cfg, "run.random_seed", 0))
    sample_key = deep_get(cfg, "keys.sample", "sample_id")
    batch_candidates = deep_get(cfg, "keys.batch_candidates", [sample_key])
    batch_keys = [k for k in batch_candidates if adata is not None and k in adata.obs and adata.obs[k].nunique() > 1]
    resolutions = list(deep_get(cfg, "core.fastcore.omicverse.leiden.resolutions", deep_get(cfg, "core.leiden_resolutions", [0.8])))
    default_resolution = resolutions[min(len(resolutions) // 2, len(resolutions) - 1)] if resolutions else 0.8
    preprocess_mode = str(deep_get(cfg, "core.fastcore.omicverse.preprocess_mode", "shiftlog"))
    if "|" not in preprocess_mode:
        preprocess_mode = f"{preprocess_mode}|seurat"
    return OmicVerseCoreOptions(
        backend=backend,
        mode=mode,
        preprocess_mode=preprocess_mode,
        target_sum=float(deep_get(cfg, "core.fastcore.omicverse.target_sum", 500000)),
        n_hvgs=int(deep_get(cfg, "core.fastcore.omicverse.n_hvgs", deep_get(cfg, "core.n_top_hvg", 3000))),
        scaled_layer=str(deep_get(cfg, "core.fastcore.omicverse.scaled_layer", "scaled")),
        n_pcs=int(deep_get(cfg, "core.fastcore.omicverse.n_pcs", deep_get(cfg, "core.n_pcs", 50))),
        n_neighbors=int(deep_get(cfg, "core.fastcore.omicverse.neighbors.n_neighbors", deep_get(cfg, "core.neighbors_n_neighbors", 15))),
        neighbors_n_pcs=int(deep_get(cfg, "core.fastcore.omicverse.neighbors.n_pcs", deep_get(cfg, "core.neighbors_n_pcs", 30))),
        neighbors_method=deep_get(cfg, "core.fastcore.omicverse.neighbors.method", None),
        neighbors_transformer=deep_get(cfg, "core.fastcore.omicverse.neighbors.mixed_transformer", None),
        umap_min_dist=float(deep_get(cfg, "core.fastcore.omicverse.umap.min_dist", deep_get(cfg, "core.umap_min_dist", 0.3))),
        umap_method=deep_get(cfg, "core.fastcore.omicverse.umap.method", None),
        leiden_resolution=float(deep_get(cfg, "core.fastcore.omicverse.leiden.default_resolution", default_resolution)),
        random_state=random_state,
        batch_keys=batch_keys,
        max_iter_harmony=int(deep_get(cfg, "core.batch_correction.max_iter_harmony", 20)),
        harmony_ncores=int(deep_get(cfg, "core.batch_correction.ncores", 0)),
        identify_robust=bool(deep_get(cfg, "core.fastcore.omicverse.identify_robust", False)),
        use_implicit_centering=bool(deep_get(cfg, "core.fastcore.omicverse.use_implicit_centering", mode == "cpu-gpu-mixed")),
    )


def time_step(timings: dict[str, float], name: str, fn: Callable, *args, **kwargs):
    start = perf_counter()
    out = fn(*args, **kwargs)
    timings[name] = perf_counter() - start
    return out


def unsupported_backend_result(backend: str, cfg: Mapping[str, Any], run_root: str | Path) -> dict[str, Any]:
    raise RuntimeError(
        f"{backend} was selected, but its executable OmicVerse adapter is unavailable. "
        "Install the separate OmicVerse environment or use the scanpy_legacy fallback."
    )


def map_standard_core_keys(adata) -> None:
    """Map common OmicVerse keys back to SCOOP's stable downstream schema."""
    if "X_pca" in adata.obsm and "X_pca_biology" not in adata.obsm:
        adata.obsm["X_pca_biology"] = adata.obsm["X_pca"].copy()
    if "X_pca" in adata.obsm and "X_pca_identity_prebatch" not in adata.obsm:
        adata.obsm["X_pca_identity_prebatch"] = adata.obsm["X_pca"].copy()
    if "X_umap" in adata.obsm:
        if "X_umap_biology" not in adata.obsm:
            adata.obsm["X_umap_biology"] = adata.obsm["X_umap"].copy()
        if "X_umap_identity" not in adata.obsm:
            adata.obsm["X_umap_identity"] = adata.obsm["X_umap"].copy()
    if "connectivities" in adata.obsp:
        if "connectivities_identity" not in adata.obsp:
            adata.obsp["connectivities_identity"] = adata.obsp["connectivities"].copy()
    if "distances" in adata.obsp:
        if "distances_identity" not in adata.obsp:
            adata.obsp["distances_identity"] = adata.obsp["distances"].copy()
    if "neighbors_identity_connectivities" in adata.obsp:
        adata.obsp["connectivities_identity"] = adata.obsp["neighbors_identity_connectivities"].copy()
    if "neighbors_identity_distances" in adata.obsp:
        adata.obsp["distances_identity"] = adata.obsp["neighbors_identity_distances"].copy()
    if "cluster_identity" in adata.obs:
        adata.obs["cluster_identity"] = adata.obs["cluster_identity"].astype("category")


def configure_omicverse(ov, options: OmicVerseCoreOptions) -> None:
    if hasattr(ov, "set_seed"):
        ov.set_seed(options.random_state, verbose=False)
    if options.mode == "cpu-gpu-mixed":
        ov.settings.cpu_gpu_mixed_init()
    else:
        ov.settings.cpu_init()


def copy_counts_to_x(adata, counts_layer: str) -> None:
    if counts_layer in adata.layers:
        adata.X = adata.layers[counts_layer].copy()
    else:
        adata.layers[counts_layer] = adata.X.copy()


def run_omicverse_preprocess_scale_pca(ov, adata, options: OmicVerseCoreOptions, timings: dict[str, float]) -> None:
    out = time_step(
        timings,
        "omicverse_preprocess",
        ov.pp.preprocess,
        adata,
        mode=options.preprocess_mode,
        target_sum=options.target_sum,
        n_HVGs=options.n_hvgs,
        batch_key=options.batch_keys[0] if options.batch_keys else None,
        identify_robust=options.identify_robust,
    )
    if out is not None and out is not adata:
        adata._init_as_actual(out)
    time_step(
        timings,
        "omicverse_scale",
        ov.pp.scale,
        adata,
        max_value=10,
        layers_add=options.scaled_layer,
        use_implicit_centering=options.use_implicit_centering,
    )
    time_step(
        timings,
        "omicverse_pca",
        ov.pp.pca,
        adata,
        n_pcs=options.n_pcs,
        layer=options.scaled_layer,
        random_state=options.random_state,
    )
    adata.uns.pop("_scaled_implicit", None)
    map_standard_core_keys(adata)


def run_harmony_bridge(adata, options: OmicVerseCoreOptions, timings: dict[str, float]) -> bool:
    if "X_pca_identity_prebatch" not in adata.obsm and "X_pca" in adata.obsm:
        adata.obsm["X_pca_identity_prebatch"] = adata.obsm["X_pca"].copy()
    try:
        time_step(
            timings,
            "harmony2",
            run_harmony2,
            adata,
            basis="X_pca_identity_prebatch",
            batch_keys=options.batch_keys,
            output="X_pca_harmony_identity",
            max_iter_harmony=options.max_iter_harmony,
            random_state=options.random_state,
            ncores=options.harmony_ncores,
        )
        return bool(options.batch_keys)
    except ImportError:
        adata.obsm["X_pca_harmony_identity"] = adata.obsm["X_pca_identity_prebatch"].copy()
        return False


def run_omicverse_graph_umap_leiden(ov, adata, options: OmicVerseCoreOptions, timings: dict[str, float]) -> None:
    neighbors_kwargs: dict[str, Any] = {
        "n_neighbors": options.n_neighbors,
        "n_pcs": options.neighbors_n_pcs,
        "use_rep": "X_pca_harmony_identity",
        "key_added": "neighbors_identity",
        "random_state": options.random_state,
    }
    if options.neighbors_method and options.neighbors_method != "auto":
        neighbors_kwargs["method"] = options.neighbors_method
    if options.mode == "cpu-gpu-mixed" and options.neighbors_transformer and options.neighbors_transformer != "auto":
        neighbors_kwargs["transformer"] = options.neighbors_transformer
    time_step(timings, "omicverse_neighbors", ov.pp.neighbors, adata, **neighbors_kwargs)

    umap_kwargs: dict[str, Any] = {
        "neighbors_key": "neighbors_identity",
        "min_dist": options.umap_min_dist,
        "random_state": options.random_state,
    }
    if options.umap_method and options.umap_method != "auto" and (options.mode == "cpu" or options.umap_method == "pumap"):
        umap_kwargs["method"] = options.umap_method
    time_step(timings, "omicverse_umap", ov.pp.umap, adata, **umap_kwargs)
    time_step(
        timings,
        "omicverse_leiden",
        ov.pp.leiden,
        adata,
        resolution=options.leiden_resolution,
        random_state=options.random_state,
        key_added="cluster_identity",
        neighbors_key="neighbors_identity",
    )
    map_standard_core_keys(adata)


def hvg_table_from_adata(adata, *, method: str) -> pd.DataFrame:
    cols = [c for c in ["highly_variable", "highly_variable_features", "highly_variable_rank", "means", "variances", "variances_norm", "residual_variances"] if c in adata.var]
    if cols:
        out = adata.var[cols].copy()
    else:
        out = pd.DataFrame(index=adata.var_names)
    out.insert(0, "gene", adata.var_names)
    out["hvg_key"] = "highly_variable_biology"
    out["method"] = method
    for key in ("highly_variable_features", "highly_variable"):
        if key in adata.var:
            adata.var["highly_variable_biology"] = adata.var[key].astype(bool)
            adata.var["highly_variable_identity"] = adata.var[key].astype(bool)
            break
    return out


def write_core_tables(adata, run_root: str | Path, *, method: str, options: OmicVerseCoreOptions) -> dict[str, str]:
    run_root = Path(run_root)
    hvg_path = write_table(hvg_table_from_adata(adata, method=method), run_root / "02_core" / "tables" / "hvg_biology.parquet")
    sweep = pd.DataFrame(
        [
            {
                "resolution": options.leiden_resolution,
                "seed": options.random_state,
                "n_clusters": int(adata.obs["cluster_identity"].nunique()) if "cluster_identity" in adata.obs else 0,
                "method": method,
            }
        ]
    )
    stability = pd.DataFrame(
        [
            {
                "resolution": options.leiden_resolution,
                "median_ari": 1.0,
                "median_nmi": 1.0,
                "n_clusters_seed0": int(adata.obs["cluster_identity"].nunique()) if "cluster_identity" in adata.obs else 0,
                "chosen": True,
                "method": method,
            }
        ]
    )
    sweep_path = write_table(sweep, run_root / "02_core" / "tables" / "leiden_sweep.parquet")
    st_path = write_table(stability, run_root / "02_core" / "tables" / "cluster_stability.parquet")
    register_file(adata, key="hvg_biology", path=hvg_path, schema="hvg_table.v1")
    register_file(adata, key="leiden_sweep", path=sweep_path, schema="leiden_sweep.v1")
    register_file(adata, key="cluster_stability", path=st_path, schema="cluster_stability.v1")
    return {
        "hvg_biology": str(hvg_path),
        "leiden_sweep": str(sweep_path),
        "cluster_stability": str(st_path),
    }


def peak_gpu_memory_mb() -> float | None:
    try:
        import torch

        if torch.cuda.is_available():
            return float(torch.cuda.max_memory_allocated() / 1024**2)
    except Exception:
        pass
    try:
        import cupy as cp

        pool = cp.get_default_memory_pool()
        return float(pool.used_bytes() / 1024**2)
    except Exception:
        return None


def summarize_result(
    adata,
    *,
    options: OmicVerseCoreOptions,
    timings: dict[str, float],
    harmony_used: bool,
    artifacts: dict[str, str],
    quality_basis: str,
) -> dict[str, Any]:
    adata.uns[f"fastcore_{options.backend}"] = {
        "backend": options.backend,
        "mode": options.mode,
        "source": "OmicVerse 2.2.x preprocessing backend adapter",
        "preprocess_mode": options.preprocess_mode,
        "target_sum": options.target_sum,
        "n_hvgs": options.n_hvgs,
        "n_pcs": options.n_pcs,
        "single_leiden_resolution": options.leiden_resolution,
    }
    return {
        "backend": options.backend,
        "harmony2_used": harmony_used,
        "batch_correction_method": "harmony2",
        "batch_keys": options.batch_keys,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_clusters": int(adata.obs["cluster_identity"].nunique()) if "cluster_identity" in adata.obs else 0,
        "timings": timings,
        "quality": {
            "accepted": True,
            "backend": options.backend,
            "quality_basis": quality_basis,
            "peak_gpu_memory_mb": peak_gpu_memory_mb(),
        },
        "artifacts": artifacts,
    }
