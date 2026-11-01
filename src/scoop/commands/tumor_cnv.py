#!/usr/bin/env python
"""Run pooled-reference FastCNVpy tumor evidence after major-lineage annotation."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path

import anndata as ad

from fastcnvpy.reference import read_gene_metadata
from scsp_agent_sop.config import deep_get, read_yaml, resolve_run_root
from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.measurement import FeatureMeasurement, feature_measurement_path
from scsp_agent_sop.storage import changed_obs_columns, ensure_dir, init_file_registry, register_file, write_obs_sidecar
from scsp_agent_sop.tumor_cnv import run_tumor_fastcnv_workflow


def main() -> dict:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--output", default=None)
    ap.add_argument("--output-format", choices=["h5ad", "obs-parquet"], default="h5ad",
                    help="obs-parquet writes only the new or changed obs columns keyed by obs_name, plus a JSON receipt")
    ap.add_argument("--gene-metadata", required=True)
    ap.add_argument("--sample-key", default=None)
    ap.add_argument("--major-lineage-key", default=None)
    ap.add_argument(
        "--parenchymal-lineages", default=None, help="comma-separated major lineage labels, e.g. Epithelial,Glial"
    )
    ap.add_argument("--normal-status-key", default=None)
    ap.add_argument("--normal-status-labels", default=None)
    ap.add_argument("--layer", default=None)
    ap.add_argument("--h5ad-mode", choices=["dense", "sparse"], default=None)
    ap.add_argument("--n-jobs", type=int, default=None)
    ap.add_argument("--reference-eligible-key")
    ap.add_argument("--compute-clusters", action="store_true", default=None, help="Opt in to quadratic CNV cell clustering")
    ap.add_argument("--cell-block-size", type=int, default=None)
    ap.add_argument("--all-cells-q75-fallback", action="store_true", default=None,
                    help="Opt in to an all-cells q75 threshold when fewer than 10 reference cells remain")
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    input_path = Path(args.input) if args.input else run_root / "artifacts" / "adata_annotation_evidence.h5ad"
    sidecar = args.output_format == "obs-parquet"
    default_output = "tumor_fastcnv_obs.parquet" if sidecar else "adata_tumor_fastcnv.h5ad"
    output_path = Path(args.output) if args.output else run_root / "artifacts" / default_output
    output_dir = run_root / "annotation" / "fastcnvpy_tumor"

    adata = ad.read_h5ad(input_path)
    init_file_registry(adata, deep_get(cfg, "run.run_id", run_root.name))
    input_obs = adata.obs.copy() if sidecar else None
    gene_metadata = read_gene_metadata(Path(args.gene_metadata))
    sample_key = args.sample_key or deep_get(cfg, "keys.sample", "sample_id")
    major_key = args.major_lineage_key or deep_get(cfg, "annotation.major_lineage_key", "cell_type_lvl1")
    parenchymal = args.parenchymal_lineages or deep_get(cfg, "tumor_fastcnv.parenchymal_lineages", "Epithelial")
    if isinstance(parenchymal, str):
        parenchymal = [item for item in parenchymal.split(",") if item]
    normal_labels = args.normal_status_labels or deep_get(cfg, "tumor_fastcnv.normal_status_labels", ["normal", "non_tumor", "healthy"])
    if isinstance(normal_labels, str):
        normal_labels = [item for item in normal_labels.split(",") if item]
    h5ad_mode = args.h5ad_mode or deep_get(cfg, "tumor_fastcnv.h5ad_mode", "sparse")
    n_jobs = args.n_jobs if args.n_jobs is not None else int(deep_get(cfg, "tumor_fastcnv.n_jobs", 1))
    cell_block_size = args.cell_block_size if args.cell_block_size is not None else int(deep_get(cfg, "tumor_fastcnv.cell_block_size", 10000))
    compute_clusters = args.compute_clusters if args.compute_clusters is not None else bool(deep_get(cfg, "tumor_fastcnv.compute_clusters", False))
    allowlist = deep_get(cfg, "tumor_fastcnv.nonparenchymal_lineages", None)
    q75_fallback = (args.all_cells_q75_fallback if args.all_cells_q75_fallback is not None
                    else bool(deep_get(cfg, "tumor_fastcnv.all_cells_q75_fallback", False)))
    measurement_path = feature_measurement_path(adata)

    result, plan = run_tumor_fastcnv_workflow(
        adata,
        gene_metadata,
        output_dir=output_dir,
        sample_key=sample_key,
        major_lineage_key=major_key,
        parenchymal_lineages=parenchymal,
        normal_status_key=args.normal_status_key or deep_get(cfg, "tumor_fastcnv.normal_status_key", None),
        normal_status_labels=normal_labels,
        layer=args.layer,
        h5ad_mode=h5ad_mode,
        n_jobs=n_jobs,
        reference_eligible_key=args.reference_eligible_key or deep_get(cfg, "tumor_fastcnv.reference_eligible_key", None),
        nonparenchymal_lineages=None if allowlist is None else set(allowlist),
        all_cells_q75_fallback=q75_fallback,
        compute_clusters=compute_clusters,
        cell_block_size=cell_block_size,
        measurement=FeatureMeasurement.load(measurement_path) if measurement_path else None,
    )

    manifest = output_dir / "tumor_fastcnv_fastcnvpy_pooled_manifest.json"
    threshold_source = str(adata.obs["fastcnv_threshold_source"].iloc[0]) if result is not None else None
    if result is not None:
        register_file(
            adata,
            key="tumor_fastcnv_manifest",
            path=manifest,
            schema="fastcnvpy.pooled_result.v1",
        )
    ensure_dir(output_path.parent)
    if sidecar:
        sidecar_receipt = write_obs_sidecar(adata, changed_obs_columns(input_obs, adata), output_path, key="tumor_fastcnv_obs")
    else:
        adata.write_h5ad(output_path)
    log_decision(
        run_root,
        module="annotation",
        decision="tumor_fastcnv_completed" if result is not None else "tumor_fastcnv_skipped",
        reason=plan.reason,
        parameters={
            "sample_key": sample_key,
            "major_lineage_key": major_key,
            "parenchymal_lineages": parenchymal,
            "reference_key": plan.reference_key,
            "reference_labels": plan.reference_labels,
            "nonparenchymal_lineages": plan.nonparenchymal_lineages,
            "all_cells_q75_fallback": q75_fallback,
            "h5ad_mode": h5ad_mode,
            "n_jobs": n_jobs,
            "cell_block_size": cell_block_size,
            "compute_clusters": compute_clusters,
        },
        evidence={
            "n_reference_cells": plan.n_reference_cells,
            "n_normal_nonparenchymal": plan.n_normal_nonparenchymal,
            "n_normal_parenchymal": plan.n_normal_parenchymal,
            "n_candidate_cells": plan.n_candidate_cells,
            "dropped_reference_by_lineage": plan.dropped_reference_by_lineage,
            "threshold_source": threshold_source,
            "measurement": plan.measurement,
            "output_obs_sidecar" if sidecar else "output_h5ad": str(output_path),
        },
    )
    limitations = []
    if result is None:
        limitations.append(f"FastCNV was not run: {plan.reason}; no CNV evidence columns were written.")
    elif threshold_source != "reference_q99_mad":
        limitations.append(f"CNV threshold is not reference-derived ({threshold_source}); "
                           "non-reference calls are unassessed or query-quantile based.")
    if plan.dropped_reference_by_lineage:
        limitations.append(f"Reviewed cells outside the nonparenchymal vocabulary were not used as reference: "
                           f"{plan.dropped_reference_by_lineage}")
    if plan.measurement["status"] == "applied" and plan.measurement["n_positioned_genes_excluded"]:
        limitations.append(f"{plan.measurement['n_positioned_genes_excluded']} positioned genes were excluded because at "
                           f"least one sample did not measure them; see {plan.measurement['excluded_genes_table']}.")
    return {"inputs": [str(input_path.resolve()), str(Path(args.gene_metadata).resolve()),
                       *([str(measurement_path.resolve())] if measurement_path else [])],
            "outputs": [str(output_path.resolve()), *([str(sidecar_receipt.resolve())] if sidecar else []),
                        *([str(manifest.resolve())] if result is not None else []),
                        *([str(Path(plan.measurement["excluded_genes_table"]).resolve())]
                          if plan.measurement["status"] == "applied" else [])],
            "execution": {"status": "completed"} if result is not None else {"status": "not_run", "reason": plan.reason},
            "validation": {"status": "not_assessed"}, "interpretation": {"status": "not_assessed"},
            "reference_plan": asdict(plan), "threshold_source": threshold_source, "output_format": args.output_format,
            "limitations": limitations}


if __name__ == "__main__":
    main()
