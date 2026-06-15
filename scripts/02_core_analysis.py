#!/usr/bin/env python
"""Run normalisation, HVG, PCA, Harmony, graph, UMAP and Leiden sweep.

Usage:
  python scripts/02_core_analysis.py --config runs/<run_id>/config/run.yaml

The script writes only final cluster and embeddings to H5AD. Sweep tables and
HVG diagnostics are written as sidecar files.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import anndata as ad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from scsp_agent_sop.config import read_yaml, deep_get, resolve_run_root
from scsp_agent_sop.core import normalize_log1p, select_hvg, run_pca, score_programs, robust_z_by_sample, build_identity_hvg_from_program_decision, run_harmony_pytorch, neighbors_umap, leiden_sweep
from scsp_agent_sop.storage import write_table, register_file, init_file_registry, ensure_dir
from scsp_agent_sop.decision_log import log_decision


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--output", default=None)
    args = ap.parse_args()
    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    input_path = Path(args.input) if args.input else run_root / "artifacts" / "adata_qc.h5ad"
    output = Path(args.output) if args.output else run_root / "artifacts" / "adata_core.h5ad"
    adata = ad.read_h5ad(input_path)
    init_file_registry(adata, deep_get(cfg, "run.run_id", run_root.name))

    normalize_log1p(adata, counts_layer=deep_get(cfg, "qc.counts_layer", "counts"), target_sum=deep_get(cfg, "core.normalize_total_target_sum", 10000))
    score_programs(adata, organism=deep_get(cfg, "run.organism", "human"), layer="log1p_norm")
    robust_z_by_sample(adata, ["stress_score", "ribo_score", "proliferation_score"], sample_key=deep_get(cfg, "keys.sample", "sample_id"))

    hvg = select_hvg(
        adata,
        counts_layer=deep_get(cfg, "qc.counts_layer", "counts"),
        batch_key=deep_get(cfg, "keys.sample", "sample_id"),
        flavor=deep_get(cfg, "core.hvg_flavor", "seurat_v3"),
        n_top_genes=deep_get(cfg, "core.n_top_hvg", 3000),
        output_key="highly_variable_biology",
    )
    hvg_path = write_table(hvg, run_root / "02_core" / "tables" / "hvg_biology.parquet")
    register_file(adata, key="hvg_biology", path=hvg_path, schema="hvg_table.v1")

    run_pca(adata, hvg_key="highly_variable_biology", obsm_key="X_pca_biology", n_comps=deep_get(cfg, "core.n_pcs", 50))
    identity_hvg = build_identity_hvg_from_program_decision(adata)
    id_path = write_table(identity_hvg, run_root / "02_core" / "tables" / "hvg_identity.parquet")
    register_file(adata, key="hvg_identity", path=id_path, schema="hvg_identity.v1")
    run_pca(adata, hvg_key="highly_variable_identity", obsm_key="X_pca_identity_prebatch", n_comps=deep_get(cfg, "core.n_pcs", 50))

    batch_keys = [k for k in deep_get(cfg, "keys.batch_candidates", ["sample_id"]) if k in adata.obs]
    try:
        run_harmony_pytorch(adata, basis="X_pca_identity_prebatch", batch_keys=batch_keys, output="X_pca_harmony_identity")
        harmony_used = True
    except ImportError:
        adata.obsm["X_pca_harmony_identity"] = adata.obsm["X_pca_identity_prebatch"].copy()
        harmony_used = False

    neighbors_umap(
        adata,
        use_rep="X_pca_harmony_identity",
        prefix="identity",
        n_neighbors=deep_get(cfg, "core.neighbors_n_neighbors", 15),
        n_pcs=None,
        min_dist=deep_get(cfg, "core.umap_min_dist", 0.3),
        random_state=deep_get(cfg, "run.random_seed", 0),
    )
    neighbors_umap(
        adata,
        use_rep="X_pca_biology",
        prefix="biology",
        n_neighbors=deep_get(cfg, "core.neighbors_n_neighbors", 15),
        n_pcs=None,
        min_dist=deep_get(cfg, "core.umap_min_dist", 0.3),
        random_state=deep_get(cfg, "run.random_seed", 0),
    )
    sweep, stability = leiden_sweep(
        adata,
        graph_prefix="identity",
        resolutions=deep_get(cfg, "core.leiden_resolutions", [0.4, 0.8, 1.2]),
        seeds=deep_get(cfg, "core.leiden_seeds", [0, 1, 2, 3, 4]),
    )
    sweep_path = write_table(sweep, run_root / "02_core" / "tables" / "leiden_sweep.parquet")
    st_path = write_table(stability, run_root / "02_core" / "tables" / "cluster_stability.parquet")
    register_file(adata, key="leiden_sweep", path=sweep_path, schema="leiden_sweep.v1")
    register_file(adata, key="cluster_stability", path=st_path, schema="cluster_stability.v1")

    ensure_dir(output.parent)
    adata.write_h5ad(output)
    log_decision(
        run_root,
        module="core_analysis",
        decision="core_analysis_complete",
        reason="Default normalization/HVG/PCA/Harmony/kNN/UMAP/Leiden sweep executed; only final clusters retained in obs.",
        parameters={"harmony_pytorch_used": harmony_used, "batch_keys": batch_keys},
        evidence={"n_clusters": int(adata.obs['cluster_identity'].nunique())},
        fallback_used=not harmony_used,
        fallback_reason=None if harmony_used else "harmony-pytorch unavailable; identity PCA used without correction",
    )


if __name__ == "__main__":
    main()
