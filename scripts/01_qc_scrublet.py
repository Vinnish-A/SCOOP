#!/usr/bin/env python
"""Run per-sample QC and Scanpy Scrublet.

Usage:
  python scripts/01_qc_scrublet.py --config runs/<run_id>/config/run.yaml

Inputs:
  - run.input_h5ad or artifacts/validated_adata.h5ad
  - layers[counts]
  - obs[sample_id]

Outputs:
  - artifacts/adata_qc.h5ad
  - 01_qc/tables/qc_thresholds_by_sample.tsv
  - 01_qc/tables/scrublet_summary_by_sample.tsv
  - 01_qc/tables/scrublet_scores.parquet
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import anndata as ad
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from scsp_agent_sop.config import read_yaml, deep_get, resolve_run_root, project_root_from_run_root
from scsp_agent_sop.qc import annotate_gene_flags, compute_basic_qc, assign_qc_flags_per_sample, ThresholdConfig, run_scrublet_per_sample
from scsp_agent_sop.storage import write_table, init_file_registry, register_file, ensure_dir
from scsp_agent_sop.decision_log import log_decision


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--output", default=None)
    args = ap.parse_args()
    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    project_root = project_root_from_run_root(run_root)
    input_path = Path(args.input) if args.input else (run_root / "artifacts" / "validated_adata.h5ad")
    if not input_path.exists():
        input_path = project_root / deep_get(cfg, "run.input_h5ad")
    output = Path(args.output) if args.output else run_root / "artifacts" / "adata_qc.h5ad"
    adata = ad.read_h5ad(input_path)
    init_file_registry(adata, deep_get(cfg, "run.run_id", run_root.name))

    organism = deep_get(cfg, "run.organism", "human")
    sample_key = deep_get(cfg, "keys.sample", "sample_id")
    counts_layer = deep_get(cfg, "qc.counts_layer", "counts")
    annotate_gene_flags(adata, organism=organism)
    compute_basic_qc(adata, counts_layer=counts_layer)

    tc = ThresholdConfig(
        low_counts_mad=deep_get(cfg, "qc.mad.low_counts", 3.0),
        low_genes_mad=deep_get(cfg, "qc.mad.low_genes", 3.0),
        high_mt_mad=deep_get(cfg, "qc.mad.high_mt", 3.0),
        high_ribo_mad=deep_get(cfg, "qc.mad.high_ribo", 4.0),
        high_hb_mad=deep_get(cfg, "qc.mad.high_hb", 4.0),
    )
    thresholds = assign_qc_flags_per_sample(adata, sample_key=sample_key, cfg=tc)
    qc_path = write_table(thresholds, run_root / "01_qc" / "tables" / "qc_thresholds_by_sample.tsv")
    register_file(adata, key="qc_thresholds_by_sample", path=qc_path, schema="qc_thresholds.v1")

    scrublet_summary = run_scrublet_per_sample(
        adata,
        sample_key=sample_key,
        counts_layer=counts_layer,
        expected_doublet_rate=deep_get(cfg, "qc.scrublet.expected_doublet_rate", 0.05),
        stdev_doublet_rate=deep_get(cfg, "qc.scrublet.stdev_doublet_rate", 0.02),
        sim_doublet_ratio=deep_get(cfg, "qc.scrublet.sim_doublet_ratio", 2.0),
        random_state=deep_get(cfg, "qc.scrublet.random_state", 0),
    )
    s_path = write_table(scrublet_summary, run_root / "01_qc" / "tables" / "scrublet_summary_by_sample.tsv")
    register_file(adata, key="scrublet_summary_by_sample", path=s_path, schema="scrublet_summary.v1")
    score_cols = [c for c in [sample_key, "doublet_score", "doublet_call_scrublet", "doublet_call"] if c in adata.obs]
    score_path = write_table(adata.obs[score_cols].reset_index(names="obs_name"), run_root / "01_qc" / "tables" / "scrublet_scores.parquet")
    register_file(adata, key="scrublet_scores", path=score_path, schema="scrublet_scores.v1")

    ensure_dir(output.parent)
    adata.write_h5ad(output)
    log_decision(
        run_root,
        module="qc",
        decision="qc_and_scrublet_complete",
        reason="QC thresholds are sample-specific; Scrublet was run per sample on raw counts.",
        parameters={"counts_layer": counts_layer, "sample_key": sample_key},
        evidence={"n_obs": int(adata.n_obs), "qc_fail_fraction": float((adata.obs['qc_class'] == 'fail').mean())},
        human_review_required=bool((thresholds["fail_fraction"] > 0.30).any()),
        review_reason="At least one sample has QC fail fraction > 30%" if (thresholds["fail_fraction"] > 0.30).any() else None,
    )


if __name__ == "__main__":
    main()
