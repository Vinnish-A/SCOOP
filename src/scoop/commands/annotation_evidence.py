#!/usr/bin/env python
"""Export cluster marker and annotation evidence tables.

Usage:
  scoop annotate evidence --config runs/<run_id>/config/run.yaml

This script is not a cell-type annotator. It exports data-derived evidence
that a subagent or analyst can use with biological knowledge before commit.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad

from scsp_agent_sop.annotation import (
    build_annotation_evidence_template,
    run_markers_fastde,
    run_markers_omicverse,
    run_markers_scanpy,
)
from scsp_agent_sop.config import deep_get, read_yaml, resolve_run_root
from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.storage import ensure_dir, init_file_registry, register_file, write_table


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--output", default=None)
    ap.add_argument("--cluster-key", default=None)
    ap.add_argument("--layer", default="log1p_norm")
    ap.add_argument("--use-fastde-cosg", action="store_true")
    ap.add_argument("--use-omicverse-cosg", action="store_true")
    args = ap.parse_args()
    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    input_path = Path(args.input) if args.input else run_root / "artifacts" / "adata_programs.h5ad"
    output = Path(args.output) if args.output else run_root / "artifacts" / "adata_annotation_evidence.h5ad"
    adata = ad.read_h5ad(input_path)
    init_file_registry(adata, deep_get(cfg, "run.run_id", run_root.name))

    groupby = args.cluster_key or deep_get(cfg, "annotation.decision.cluster_key", "cluster_identity")
    if args.use_fastde_cosg:
        markers = run_markers_fastde(
            adata,
            groupby=groupby,
            layer=args.layer,
            method="cosg",
            n_genes=deep_get(cfg, "annotation.marker_n_genes", 100),
        )
        method = "fastde_cosg"
    elif args.use_omicverse_cosg:
        markers = run_markers_omicverse(
            adata,
            groupby=groupby,
            layer=args.layer,
            method="cosg",
            n_genes=deep_get(cfg, "annotation.marker_n_genes", 100),
        )
        method = "omicverse_cosg"
    else:
        markers = run_markers_scanpy(
            adata,
            groupby=groupby,
            layer=args.layer,
            method=deep_get(cfg, "annotation.default_marker_method", "wilcoxon"),
            n_genes=deep_get(cfg, "annotation.marker_n_genes", 100),
        )
        method = "scanpy"
    marker_path = write_table(markers, run_root / "annotation" / "tables" / "cluster_markers.parquet")
    group_column = next((name for name in ("group", "cluster", "cell_type") if name in markers.columns), None)
    gene_column = next((name for name in ("names", "gene", "feature") if name in markers.columns), None)
    if group_column and gene_column:
        for cluster_id, frame in markers.groupby(group_column, observed=True):
            adata.uns[f"top_markers_{cluster_id}"] = frame[gene_column].dropna().astype(str).head(100).tolist()
    evidence = build_annotation_evidence_template(adata, cluster_key=groupby)
    evidence_path = write_table(evidence, run_root / "annotation" / "tables" / "annotation_evidence_template.tsv")
    register_file(adata, key="cluster_markers", path=marker_path, schema="cluster_markers.v1")
    register_file(adata, key="annotation_evidence_template", path=evidence_path, schema="annotation_evidence.v1")

    ensure_dir(output.parent)
    adata.write_h5ad(output)
    log_decision(
        run_root,
        module="annotation",
        decision="annotation_evidence_exported",
        reason="Marker table and evidence template exported; no biological labels were assigned.",
        parameters={"marker_method": method},
    )


if __name__ == "__main__":
    main()
