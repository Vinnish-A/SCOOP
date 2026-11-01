#!/usr/bin/env python
"""Run per-sample QC and Scanpy Scrublet.

Usage:
  scoop qc --config runs/<run_id>/config/run.yaml

Inputs:
  - run.input_h5ad or artifacts/validated_adata.h5ad
  - layers[counts]
  - obs[sample_id]

QC-owned obs columns already in the input (qc_*, doublet_*, pct_counts_*,
log10_*, total_counts, final_use, ...) are renamed to prior_* first; the new
final_use is ANDed with prior_final_use, so QC never re-admits a cell.

Outputs:
  - artifacts/adata_qc.h5ad, or with --output-format obs-parquet only the QC obs
    columns (with any prior_* they displaced) at --output plus a <stem>.json receipt (var gene flags and
    uns are not written)
  - qc/tables/qc_thresholds_by_sample.tsv
  - qc/tables/scrublet_summary_by_sample.tsv
  - qc/tables/scrublet_scores.parquet
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad

from scsp_agent_sop.assays import assay_profile, spatial_qc
from scsp_agent_sop.config import deep_get, project_root_from_run_root, read_yaml, resolve_obs_keys, resolve_run_root
from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.measurement import FeatureMeasurement, feature_measurement_path
from scsp_agent_sop.qc import (
    ThresholdConfig,
    annotate_gene_flags,
    assign_qc_flags_per_sample,
    compute_basic_qc,
    inherit_qc,
    is_qc_owned,
    retire_prior_qc,
    run_scrublet_per_sample,
)
from scsp_agent_sop.storage import (write_h5ad_atomic, ensure_dir, init_file_registry, register_file, write_obs_sidecar,
                                    write_table)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--output", default=None)
    ap.add_argument("--inherit-from", help="Parent QC H5AD; verify exact retained-cell counts and inherit decisions")
    ap.add_argument("--output-format", choices=["h5ad", "obs-parquet"], default="h5ad",
                    help="obs-parquet writes only the QC obs columns (and any prior_* they displaced) keyed by obs_name, "
                         "plus a JSON receipt")
    args = ap.parse_args()
    if args.output_format == "obs-parquet" and not args.output:
        # the managed stage receipt checks the --output path
        ap.error("--output-format obs-parquet needs an explicit --output <path>.parquet")
    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    project_root = project_root_from_run_root(run_root)
    input_path = Path(args.input) if args.input else (run_root / "artifacts" / "validated_adata.h5ad")
    if not input_path.exists():
        input_path = project_root / deep_get(cfg, "run.input_h5ad")
    output = Path(args.output) if args.output else run_root / "artifacts" / "adata_qc.h5ad"
    adata = ad.read_h5ad(input_path)
    init_file_registry(adata, deep_get(cfg, "run.run_id", run_root.name))

    def write_output() -> None:
        ensure_dir(output.parent)
        if args.output_format == "h5ad":
            write_h5ad_atomic(adata, output)
        else:
            # The QC-owned columns, plus the prior_* columns they displaced.
            columns = [c for c in adata.obs.columns if is_qc_owned(c.removeprefix("prior_"))]
            write_obs_sidecar(adata, columns, output, base=input_path, key="qc_obs")

    organism = deep_get(cfg, "run.organism", "human")
    obs_keys = resolve_obs_keys(cfg, adata.obs, required=("sample",))
    sample_key = obs_keys.sample
    counts_layer = deep_get(cfg, "qc.counts_layer", "counts")
    if args.inherit_from:
        provenance = inherit_qc(adata, Path(args.inherit_from), counts_layer=counts_layer, sample_key=sample_key)
        write_output()
        log_decision(run_root, module="qc", decision="qc_inherited", reason="Exact retained subset; no QC or doublet refit.",
                     parameters=provenance, evidence={"output_h5ad": str(output)})
        return
    annotate_gene_flags(adata, organism=organism, symbol_key=deep_get(cfg, "qc.symbol_key"))
    # S026: keep what earlier steps decided instead of overwriting it.
    prior = retire_prior_qc(adata)

    def keep_prior_exclusions() -> int:
        if "final_use" not in prior:
            return 0
        kept = adata.obs["prior_final_use"].to_numpy(bool)
        adata.obs["final_use"] = adata.obs["final_use"].to_numpy(bool) & kept
        return int((~kept).sum())

    prior_record = {"prior_qc_columns": prior}
    assay = assay_profile(cfg)
    if assay != "scrna":
        thresholds = spatial_qc(adata, cfg)
        prior_record["prior_excluded_cells"] = keep_prior_exclusions()
        qc_path = write_table(thresholds, run_root / "qc" / "tables" / "qc_thresholds_by_sample.tsv")
        register_file(adata, key="qc_thresholds_by_sample", path=qc_path, schema="spatial_qc_thresholds.v1")
        fov_key = deep_get(cfg, "qc.spatial.fov_key", "fov")
        if fov_key in adata.obs:
            by_fov = (
                adata.obs.groupby([sample_key, fov_key], observed=True)
                .agg(
                    n_obs=("final_use", "size"),
                    retained=("final_use", "sum"),
                    median_counts=("total_counts", "median"),
                    median_genes=("n_genes_by_counts", "median"),
                    median_non_gene_fraction=("non_gene_feature_fraction", "median"),
                )
                .reset_index()
            )
            fov_path = write_table(by_fov, run_root / "qc/tables/spatial_qc_by_fov.tsv")
            register_file(adata, key="spatial_qc_by_fov", path=fov_path, schema="spatial_fov_qc.v1")
        write_output()
        log_decision(
            run_root,
            module="qc",
            decision="spatial_qc_complete",
            reason="Explicit spatial thresholds; no droplet Scrublet or mitochondrial exclusion.",
            parameters={"assay": assay, "counts_layer": counts_layer, "sample_key": sample_key,
                        "obs_keys": obs_keys.record()},
            evidence={"n_obs": int(adata.n_obs), "retained": int(adata.obs.final_use.sum()), **prior_record},
            human_review_required=bool((thresholds.excluded_fraction > 0.30).any()),
            review_reason="Spatial exclusion exceeds 30%; inspect tissue/FOV and cell-type bias."
            if (thresholds.excluded_fraction > 0.30).any()
            else None,
        )
        return
    measurement_path = feature_measurement_path(adata)
    measurement = FeatureMeasurement.load(measurement_path) if measurement_path else None
    if measurement is not None:
        measurement.check(adata)
    compute_basic_qc(adata, counts_layer=counts_layer, measurement=measurement)

    tc = ThresholdConfig(
        low_counts_mad=deep_get(cfg, "qc.mad.low_counts", 3.0),
        low_genes_mad=deep_get(cfg, "qc.mad.low_genes", 3.0),
        high_mt_mad=deep_get(cfg, "qc.mad.high_mt", 3.0),
        high_ribo_mad=deep_get(cfg, "qc.mad.high_ribo", 4.0),
        high_hb_mad=deep_get(cfg, "qc.mad.high_hb", 4.0),
    )
    thresholds = assign_qc_flags_per_sample(adata, sample_key=sample_key, cfg=tc)
    qc_path = write_table(thresholds, run_root / "qc" / "tables" / "qc_thresholds_by_sample.tsv")
    register_file(adata, key="qc_thresholds_by_sample", path=qc_path, schema="qc_thresholds.v1")

    scrublet_summary = run_scrublet_per_sample(
        adata,
        sample_key=sample_key,
        counts_layer=counts_layer,
        expected_doublet_rate=deep_get(cfg, "qc.scrublet.expected_doublet_rate", 0.05),
        stdev_doublet_rate=deep_get(cfg, "qc.scrublet.stdev_doublet_rate", 0.02),
        sim_doublet_ratio=deep_get(cfg, "qc.scrublet.sim_doublet_ratio", 2.0),
        random_state=deep_get(cfg, "qc.scrublet.random_state", 0),
        workers=deep_get(cfg, "qc.scrublet.workers", "auto"),
        threads_per_worker=int(deep_get(cfg, "qc.scrublet.threads_per_worker", 1)),
        checkpoint_dir=run_root / "qc" / "scrublet_checkpoints",
    )
    s_path = write_table(scrublet_summary, run_root / "qc" / "tables" / "scrublet_summary_by_sample.tsv")
    register_file(adata, key="scrublet_summary_by_sample", path=s_path, schema="scrublet_summary.v1")
    score_cols = [c for c in [sample_key, "doublet_score", "doublet_call_scrublet", "doublet_call"] if c in adata.obs]
    score_path = write_table(
        adata.obs[score_cols].reset_index(names="obs_name"), run_root / "qc" / "tables" / "scrublet_scores.parquet"
    )
    register_file(adata, key="scrublet_scores", path=score_path, schema="scrublet_scores.v1")

    prior_record["prior_excluded_cells"] = keep_prior_exclusions()
    write_output()
    log_decision(
        run_root,
        module="qc",
        decision="qc_and_scrublet_complete",
        reason="QC thresholds are sample-specific; Scrublet was run per sample on raw counts.",
        parameters={
            "counts_layer": counts_layer,
            "sample_key": sample_key,
            "obs_keys": obs_keys.record(),
            "symbol_key": deep_get(cfg, "qc.symbol_key"),
            "measurement": str(measurement_path) if measurement_path else "not_declared",
        },
        evidence={
            "n_obs": int(adata.n_obs),
            "qc_fail_fraction": float((adata.obs["qc_class"] == "fail").mean()),
            "unassessed_cells": {m: int(adata.obs[f"{m}_unassessed"].sum()) for m in ("mt", "ribo", "hb")},
            "retained": int(adata.obs.final_use.sum()),
            **prior_record,
        },
        human_review_required=bool((thresholds["fail_fraction"] > 0.30).any()),
        review_reason="At least one sample has QC fail fraction > 30%"
        if (thresholds["fail_fraction"] > 0.30).any()
        else None,
    )


if __name__ == "__main__":
    main()
