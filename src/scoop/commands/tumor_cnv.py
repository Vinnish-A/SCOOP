#!/usr/bin/env python
"""Run pooled-reference FastCNVpy tumor evidence after major-lineage annotation."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad

from fastcnvpy.reference import read_gene_metadata
from scsp_agent_sop.config import deep_get, read_yaml, resolve_run_root
from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.storage import ensure_dir, init_file_registry, register_file
from scsp_agent_sop.tumor_cnv import run_tumor_fastcnv_workflow


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--output", default=None)
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
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    input_path = Path(args.input) if args.input else run_root / "artifacts" / "adata_annotation_evidence.h5ad"
    output_path = Path(args.output) if args.output else run_root / "artifacts" / "adata_tumor_fastcnv.h5ad"
    output_dir = run_root / "annotation" / "fastcnvpy_tumor"

    adata = ad.read_h5ad(input_path)
    init_file_registry(adata, deep_get(cfg, "run.run_id", run_root.name))
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
        compute_clusters=compute_clusters,
        cell_block_size=cell_block_size,
    )

    if result is not None:
        register_file(
            adata,
            key="tumor_fastcnv_manifest",
            path=output_dir / "tumor_fastcnv_fastcnvpy_pooled_manifest.json",
            schema="fastcnvpy.pooled_result.v1",
        )
    ensure_dir(output_path.parent)
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
            "output_h5ad": str(output_path),
        },
    )


if __name__ == "__main__":
    main()
