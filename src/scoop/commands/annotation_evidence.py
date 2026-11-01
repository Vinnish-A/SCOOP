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
import pandas as pd

from scsp_agent_sop.annotation import (
    build_annotation_evidence_template,
    run_markers_fastde,
    run_markers_omicverse,
    run_markers_scanpy,
)
from scsp_agent_sop.config import deep_get, read_yaml, resolve_run_root
from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.measurement import FORMAT as MEASUREMENT_FORMAT, FeatureMeasurement, feature_measurement_path
from scsp_agent_sop.storage import ensure_dir, init_file_registry, register_file, sha256_file, write_table, write_h5ad_atomic


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--output", default=None)
    ap.add_argument("--cluster-key", default=None)
    ap.add_argument("--layer", default="log1p_norm")
    ap.add_argument("--use-fastde-cosg", action="store_true")
    ap.add_argument("--use-omicverse-cosg", action="store_true")
    ap.add_argument("--feature-measurement", default=None,
                    help="scoop.feature-measurement.v1 .npz; default: uns file_registry artifacts.feature_measurement")
    args = ap.parse_args()
    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    input_path = Path(args.input) if args.input else run_root / "artifacts" / "adata_programs.h5ad"
    output = Path(args.output) if args.output else run_root / "artifacts" / "adata_annotation_evidence.h5ad"
    adata = ad.read_h5ad(input_path)
    init_file_registry(adata, deep_get(cfg, "run.run_id", run_root.name))

    decision_key = deep_get(cfg, "annotation.decision.cluster_key", "cluster_identity")
    groupby = args.cluster_key or decision_key
    # One DEG table per partition: a subclustering must not overwrite its parent's evidence.
    marker_key = "cluster_markers" if groupby == decision_key else f"cluster_markers__{groupby}"
    method = deep_get(cfg, "annotation.default_marker_method", "wilcoxon")
    if method not in {"wilcoxon", "t-test", "t-test_overestim_var"}:
        raise ValueError("annotation.default_marker_method must supply gene-level statistical tests; COSG is an optional specificity table")
    measurement_path = feature_measurement_path(adata, args.feature_measurement)
    measurement = FeatureMeasurement.load(measurement_path) if measurement_path else None
    markers = run_markers_scanpy(adata, groupby=groupby, layer=args.layer, method=method, n_genes=adata.n_vars,
                                 measurement=measurement)
    measurement_record = {"status": "not_declared"}
    if measurement is not None:
        if args.feature_measurement:
            register_file(adata, key="feature_measurement", path=measurement_path, category="artifacts", schema=MEASUREMENT_FORMAT)
        n_full = int(measurement.measured_in_all(adata.obs[measurement.sample_key]).sum())
        measurement_record = {"status": "applied", "path": str(measurement_path.resolve()), "sha256": sha256_file(measurement_path),
                              "sample_key": measurement.sample_key, "n_genes_measured_in_all_samples": n_full,
                              "n_genes_masked": adata.n_vars - n_full}
    extra_outputs = []
    if args.use_fastde_cosg or args.use_omicverse_cosg:
        runner = run_markers_fastde if args.use_fastde_cosg else run_markers_omicverse
        specificity = runner(
            adata,
            groupby=groupby,
            layer=args.layer,
            method="cosg",
            n_genes=deep_get(cfg, "annotation.marker_n_genes", 100),
            measurement=measurement,
        )
        path = write_table(specificity, run_root / "annotation/tables/cluster_specificity_scores.parquet")
        scored = {"gene_axis": "all_input_genes"}
        if measurement is not None and measurement_record["n_genes_masked"]:
            excluded = adata.var_names[~measurement.measured_in_all(adata.obs[measurement.sample_key])]
            excluded_path = write_table(pd.DataFrame({"gene": excluded}),
                                        run_root / "annotation/tables/cluster_specificity_excluded_genes.tsv")
            extra_outputs.append(str(excluded_path.resolve()))
            scored = {"gene_axis": "genes_measured_in_every_sample", "n_genes_excluded": len(excluded),
                      "excluded_genes_table": str(excluded_path.resolve())}
        register_file(adata, key="cluster_specificity_scores", path=path, schema="cluster_specificity_scores.v1", extra=scored)
        extra_outputs.append(str(path.resolve()))
        measurement_record["specificity"] = scored
    # The evidence bundle reads per-unit markers from this label-keyed table.
    marker_path = write_table(markers, run_root / "annotation" / "tables" / f"{marker_key}.parquet")
    evidence = build_annotation_evidence_template(adata, cluster_key=groupby)
    evidence_path = write_table(evidence, run_root / "annotation" / "tables" / "annotation_evidence_template.tsv")
    register_file(adata, key=marker_key, path=marker_path, schema="cluster_markers.v1", extra={"label_key": groupby})
    register_file(adata, key="annotation_evidence_template", path=evidence_path, schema="annotation_evidence.v1")

    ensure_dir(output.parent)
    write_h5ad_atomic(adata, output)
    log_decision(
        run_root,
        module="annotation",
        decision="annotation_evidence_exported",
        reason="Marker table and evidence template exported; no biological labels were assigned.",
        parameters={"marker_method": method, "measurement": measurement_record["status"]},
    )
    return {"inputs": [str(input_path.resolve()), *([measurement_record["path"]] if measurement is not None else [])],
            "outputs": [str(output.resolve()), str(marker_path.resolve()), str(evidence_path.resolve()), *extra_outputs],
            "execution": {"status": "completed"},
            "validation": {"status": "not_assessed", "scope": "exported_gene_evidence_not_biological_adjudication"},
            "interpretation": {"status": "not_assessed"},
            "gene_evidence": {"method": method, "cluster_key": groupby, "layer": args.layer,
                              "requested_genes_per_group": adata.n_vars, "exported_gene_rows": len(markers)},
            "measurement": measurement_record,
            "coverage": [{"scope": "cluster_gene_evidence", "unit": "observation_views",
                          "total": adata.n_obs, "processed": adata.n_obs}]}


if __name__ == "__main__":
    main()
