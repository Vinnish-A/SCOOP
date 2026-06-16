#!/usr/bin/env python
"""Export cluster marker and annotation evidence tables.

Usage:
  python scripts/04a_annotation_evidence.py --config runs/<run_id>/config/run.yaml

This script is not a cell-type annotator. It exports data-derived evidence
that a subagent or analyst can use with biological knowledge before commit.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import anndata as ad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from scsp_agent_sop.config import read_yaml, deep_get, resolve_run_root
from scsp_agent_sop.annotation import run_markers_scanpy, run_markers_omicverse, run_markers_fastde, build_annotation_evidence_template
from scsp_agent_sop.storage import write_table, register_file, init_file_registry, ensure_dir
from scsp_agent_sop.decision_log import log_decision


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--output", default=None)
    ap.add_argument("--use-fastde-cosg", action="store_true")
    ap.add_argument("--use-omicverse-cosg", action="store_true")
    args = ap.parse_args()
    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    input_path = Path(args.input) if args.input else run_root / "artifacts" / "adata_programs.h5ad"
    output = Path(args.output) if args.output else run_root / "artifacts" / "adata_annotation_evidence.h5ad"
    adata = ad.read_h5ad(input_path)
    init_file_registry(adata, deep_get(cfg, "run.run_id", run_root.name))

    groupby = "cluster_identity"
    if args.use_fastde_cosg:
        markers = run_markers_fastde(adata, groupby=groupby, method="cosg", n_genes=deep_get(cfg, "annotation.marker_n_genes", 100))
        method = "fastde_cosg"
    elif args.use_omicverse_cosg:
        markers = run_markers_omicverse(adata, groupby=groupby, method="cosg", n_genes=deep_get(cfg, "annotation.marker_n_genes", 100))
        method = "omicverse_cosg"
    else:
        markers = run_markers_scanpy(adata, groupby=groupby, method=deep_get(cfg, "annotation.default_marker_method", "wilcoxon"), n_genes=deep_get(cfg, "annotation.marker_n_genes", 100))
        method = "scanpy"
    marker_path = write_table(markers, run_root / "04_annotation" / "tables" / "cluster_markers.parquet")
    evidence = build_annotation_evidence_template(adata, cluster_key=groupby)
    evidence_path = write_table(evidence, run_root / "04_annotation" / "tables" / "annotation_evidence_template.tsv")
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
