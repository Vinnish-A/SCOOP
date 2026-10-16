"""Explicit spatial expression baseline, independent of the scRNA recipe."""

from __future__ import annotations

from pathlib import Path
from time import perf_counter
from typing import Any, Mapping

import numpy as np
import pandas as pd
from scipy import sparse

from scsp_agent_sop.assays import assay_profile, gene_expression_mask
from scsp_agent_sop.config import deep_get
from scsp_agent_sop.core import leiden_search_options, leiden_sweep, run_harmony2, select_hvg
from scsp_agent_sop.storage import register_file, write_table


def run_spatial_core(adata, cfg: Mapping[str, Any], run_root: str | Path) -> dict[str, Any]:
    """Compute a declared baseline; it does not assign biological labels/domains."""
    import scanpy as sc

    options = deep_get(cfg, "core.spatial", {})
    feature_method = options.get("features", "all_expressed")
    normalization = options.get("normalization", "library_log1p")
    if feature_method not in {"all_expressed", "seurat"}:
        raise ValueError(
            "Set core.spatial.features to all_expressed or seurat; "
            "panel size alone does not establish an optimal method."
        )
    if normalization not in {"library_log1p", "log1p"}:
        raise ValueError(
            "Set core.spatial.normalization to library_log1p or log1p for the intended expression/domain analysis."
        )
    layer = deep_get(cfg, "qc.counts_layer", "counts")
    counts = adata.layers[layer]
    values = counts.data if sparse.issparse(counts) else np.asarray(counts)
    if (
        not np.isfinite(values).all()
        or (values < 0).any()
        or not np.allclose(values, np.rint(values), rtol=0, atol=1e-6)
    ):
        raise ValueError("Spatial core requires finite nonnegative integer counts.")
    genes = gene_expression_mask(adata, cfg)
    eligible = genes & (np.asarray((counts > 0).sum(axis=0)).ravel() >= int(options.get("min_cells", 3)))
    if eligible.sum() < 2 or adata.n_obs < 3:
        raise ValueError("Spatial embedding needs at least three observations and two expressed target genes.")
    started = perf_counter()
    work = adata[:, genes].copy()
    work.X = work.layers[layer].copy()
    if normalization == "library_log1p":
        sc.pp.normalize_total(work, target_sum=float(options.get("target_sum", 10000)))
    work.uns.pop("log1p", None)
    sc.pp.log1p(work)
    adata.uns["log1p"] = dict(work.uns["log1p"])
    # Controls stay in the original counts. Their normalized values are zero
    # and they can never enter PCA or the feature-selection denominator.
    if sparse.issparse(counts):
        full = sparse.csr_matrix((adata.n_obs, adata.n_vars), dtype=np.float32)
        coo = sparse.coo_matrix(work.X)
        full = sparse.csr_matrix((coo.data, (coo.row, np.flatnonzero(genes)[coo.col])), shape=full.shape)
    else:
        full = np.zeros(adata.shape, dtype=np.float32)
        full[:, genes] = work.X
    adata.layers["log1p_norm"] = full
    adata.X = full.copy()
    work.layers["log1p_norm"] = work.X.copy()
    if feature_method == "seurat":
        # Dispersion-HVG works on explicitly log-transformed expression. Never
        # use count-based V3/LOESS by accident on a targeted panel.
        table = select_hvg(
            work,
            counts_layer=layer,
            batch_key=options.get("hvg_batch_key"),
            flavor="seurat",
            n_top_genes=int(options.get("n_top_genes", 2000)),
            min_cells=int(options.get("min_cells", 3)),
        )
        selected = np.zeros(adata.n_vars, bool)
        selected[genes] = work.var.highly_variable_biology.to_numpy(bool)
        adata.uns["hvg_selection"] = work.uns["hvg_selection"]
        table = table.set_index("gene").reindex(adata.var_names).reset_index()
    else:
        selected = eligible
        table = pd.DataFrame({"gene": adata.var_names, "method": "all_expressed", "highly_variable": selected})
    for key in ("highly_variable_biology", "highly_variable_identity"):
        adata.var[key] = selected
    table["selected_for_pca"] = selected
    feature_path = write_table(table, Path(run_root) / "core/tables/hvg_biology.parquet")
    register_file(adata, key="hvg_biology", path=feature_path, schema="spatial_features.v1")
    # Subset only for PCA; all native genes/counts remain in the delivered object.
    work = adata[:, selected].copy()
    work.X = work.layers["log1p_norm"].copy()
    if options.get("scale", True):
        sc.pp.scale(work, zero_center=True, max_value=10)
    # Source objects may carry a stale generic highly_variable mask.
    sc.pp.pca(
        work,
        n_comps=min(int(deep_get(cfg, "core.n_pcs", 40)), work.n_obs - 1, work.n_vars - 1),
        mask_var=None,
        random_state=int(deep_get(cfg, "run.random_seed", 0)),
    )
    for key in ("X_pca_biology", "X_pca_identity_prebatch"):
        adata.obsm[key] = work.obsm["X_pca"].copy()
    batch_method = options.get("batch_correction", "none")
    if batch_method == "harmony2":
        run_harmony2(adata, batch_keys=deep_get(cfg, "keys.batch_candidates", []))
    elif batch_method == "none":
        adata.obsm["X_pca_harmony_identity"] = adata.obsm["X_pca_identity_prebatch"].copy()
    else:
        raise ValueError("core.spatial.batch_correction must be none or harmony2.")
    seed = int(deep_get(cfg, "run.random_seed", 0))
    sc.pp.neighbors(
        adata,
        use_rep="X_pca_harmony_identity",
        key_added="neighbors_identity",
        n_neighbors=min(int(deep_get(cfg, "core.neighbors_n_neighbors", 15)), adata.n_obs - 1),
        random_state=seed,
    )
    for key in ("connectivities", "distances"):
        adata.obsp[f"{key}_identity"] = adata.obsp[f"neighbors_identity_{key}"].copy()
        adata.obsp[f"{key}_biology"] = adata.obsp[f"{key}_identity"].copy()
    sc.tl.umap(
        adata,
        neighbors_key="neighbors_identity",
        random_state=seed,
        min_dist=float(deep_get(cfg, "core.umap_min_dist", 0.3)),
    )
    for key in ("X_umap_identity", "X_umap_biology"):
        adata.obsm[key] = adata.obsm["X_umap"].copy()
    # Retain graph before clustering so a later failure does not repeat it.
    checkpoint = Path(run_root) / "artifacts/adata_graph.h5ad"
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(checkpoint)
    search = leiden_search_options(cfg, run_root)
    # No change of resolution grid, seeds, or stability strategy for spatial data.
    sweep, stability = leiden_sweep(
        adata,
        resolutions=deep_get(cfg, "core.leiden_resolutions", [0.25, 0.5, 0.75, 1.0, 1.25, 1.5]),
        seeds=deep_get(cfg, "core.leiden_seeds", [0, 1, 2, 3, 4]),
        search_config=search,
    )
    artifacts = {"hvg_biology": str(feature_path), "graph_checkpoint": str(checkpoint)}
    for key, frame in (("leiden_sweep", sweep), ("cluster_stability", stability)):
        path = write_table(frame, Path(run_root) / f"core/tables/{key}.parquet")
        register_file(adata, key=key, path=path, schema=f"{key}.v1")
        artifacts[key] = str(path)
    adata.uns["spatial_preprocessing"] = {
        "assay": assay_profile(cfg),
        "features": feature_method,
        "normalization": normalization,
        "target_sum": float(options.get("target_sum", 10000)) if normalization == "library_log1p" else 0.0,
        "scale": bool(options.get("scale", True)),
        "batch_correction": batch_method,
        "n_pca_features": int(selected.sum()),
        "annotation_status": "unreviewed",
    }
    return {
        "backend": "spatial_scanpy",
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
        "n_clusters": int(adata.obs.cluster_identity.nunique()),
        "artifacts": artifacts,
        "timings": {"spatial_total": perf_counter() - started},
        "batch_correction_method": batch_method,
        "quality": {"execution_complete": True, "biological_review": "not_assessed"},
    }
