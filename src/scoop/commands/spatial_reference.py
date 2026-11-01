#!/usr/bin/env python
"""Build a measurement-aware RCTD reference package (scoop.reference-profiles.v2).

Usage:
  scoop spatial build-reference --h5ad atlas.h5ad --label-key cell_type_lvl3 --donor-key donor_id --output ref/

Counts are streamed; only reviewed, non-technical cells enter (see --status-key / --admit-status).
Each component's review status is the most cautious of SCOOP's committed {broad,fine}_annotation_status
at --annotation-level and --review-status-key; with neither it is evidence_pending, which production
fits refuse. A FeatureMeasurement (--feature-measurement, else the one registered in the H5AD) gives per-gene
measured-cell denominators. Pass the package to `scoop spatial` through spatial.rctd.reference_profiles
and reference_manifest; its heldout_validation_genes.tsv lists the genes for spatial.rctd.exclude_genes.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from scsp_agent_sop.measurement import FeatureMeasurement
from scsp_agent_sop.reference_profiles import build_reference, registered_measurement


def _table(path: str) -> pd.DataFrame:
    table = pd.read_parquet(path) if path.endswith(".parquet") else pd.read_csv(path, sep="\t")
    return table.set_index("obs_name") if "obs_name" in table else table


def main() -> dict:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", required=True, help="native integer counts (X or --layer)")
    ap.add_argument("--output", required=True, help="package directory")
    ap.add_argument("--label-key", required=True, help="obs column naming each cell's component")
    ap.add_argument("--donor-key", required=True)
    ap.add_argument("--labels", help="parquet/TSV joined to obs by obs name (index or obs_name column)")
    ap.add_argument("--layer", default=None)
    ap.add_argument("--status-key", default="annotation_status")
    ap.add_argument("--admit-status", nargs="+", default=["accepted", "labeled"])
    ap.add_argument("--reviewed-key", default="annotation_reviewed", help="must be true when the column exists")
    ap.add_argument("--technical-key", help="extra boolean column of technical cells, never admitted")
    ap.add_argument("--require-true", nargs="*", default=[], help="boolean obs columns that must be true (e.g. final_use)")
    ap.add_argument("--review-status-key", help="extra per-cell review status column; can only add caution to the committed status")
    ap.add_argument("--annotation-level", default="fine", choices=["fine", "broad"],
                    help="annotation level of --label-key; selects the committed {level}_annotation_status")
    ap.add_argument("--parent-key")
    ap.add_argument("--compartment-key")
    ap.add_argument("--components", help="TSV: component[, parent, compartment, review_status]; overrides per-cell "
                                         "parent/compartment; review_status can only add caution")
    ap.add_argument("--group-key", help="keep linear sums per group (e.g. dataset) for exact re-aggregation")
    ap.add_argument("--chemistry-key", help="obs column for mean native exposure by chemistry")
    ap.add_argument("--exposure-key", help="obs column of full native exposure; default: gene-feature row sums (extra "
                                           "pass), checked against obs native_full_gene_UMI when present")
    ap.add_argument("--feature-measurement", help="scoop.feature-measurement.v1 .npz")
    ap.add_argument("--min-umi", type=int, default=100)
    ap.add_argument("--min-cells", type=int, default=25)
    ap.add_argument("--min-donors", type=int, default=2)
    ap.add_argument("--min-measured-fraction", type=float, default=0.5)
    ap.add_argument("--min-measured-cells", type=int, default=25)
    ap.add_argument("--holdout-fraction", type=float, default=0.1)
    ap.add_argument("--holdout-max-cells", type=int, default=200)
    ap.add_argument("--holdout-min-cells", type=int, default=20)
    ap.add_argument("--heldout-genes-per-component", type=int, default=3)
    ap.add_argument("--heldout-random-genes", type=int, default=50)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--chunk-rows", type=int, default=50000)
    ap.add_argument("--seed", type=int, default=20261030)
    ap.add_argument("--config", help="run config; only places the stage receipt under its run root")
    args = ap.parse_args()
    measurement_path = Path(args.feature_measurement) if args.feature_measurement else registered_measurement(args.h5ad)
    manifest = build_reference(
        args.h5ad, args.output, label_key=args.label_key, donor_key=args.donor_key,
        labels=_table(args.labels) if args.labels else None, layer=args.layer, status_key=args.status_key,
        admit_status=args.admit_status, reviewed_key=args.reviewed_key, technical_key=args.technical_key,
        require_true=args.require_true, review_status_key=args.review_status_key, annotation_level=args.annotation_level,
        parent_key=args.parent_key,
        compartment_key=args.compartment_key, components=pd.read_csv(args.components, sep="\t") if args.components else None,
        group_key=args.group_key, chemistry_key=args.chemistry_key, exposure_key=args.exposure_key,
        measurement=FeatureMeasurement.load(measurement_path) if measurement_path else None,
        min_umi=args.min_umi, min_cells=args.min_cells, min_donors=args.min_donors,
        min_measured_fraction=args.min_measured_fraction, min_measured_cells=args.min_measured_cells,
        holdout_fraction=args.holdout_fraction, holdout_max_cells=args.holdout_max_cells,
        holdout_min_cells=args.holdout_min_cells, heldout_genes_per_component=args.heldout_genes_per_component,
        heldout_random_genes=args.heldout_random_genes, workers=args.workers, chunk_rows=args.chunk_rows, seed=args.seed)
    out = Path(args.output).resolve()
    holdout_ok = all(c["holdout_status"] in {"sufficient", "not_requested"} for c in manifest["components"])
    print(out / "reference_manifest.json")
    # The counts H5AD is bound in the manifest by path, size and var hash; hashing an atlas-scale
    # input here would cost another full read.
    return {"inputs": [p for p in [args.labels, args.components, measurement_path and str(measurement_path)] if p],
            "outputs": [str(out / name) for name in ["reference_manifest.json", "reference_profiles.npz", "component_sums.npz",
                                                     "holdout_cells.h5ad", "heldout_validation_genes.tsv", "reference_coverage.tsv"]],
            "execution": {"status": "completed"},
            "validation": {"status": "limited", "scope": "integer_counts_admission_and_measurement_denominators",
                           "checks": [{"name": "integer_native_counts", "status": "passed"},
                                      {"name": "measured_cell_denominators",
                                       "status": "passed" if measurement_path else "not_applicable"},
                                      {"name": "holdout_rule_met_by_every_component", "status": "passed" if holdout_ok else "limited"}]},
            "interpretation": {"status": "not_applicable"},
            "coverage": [{"scope": "admitted_reference_cells", "unit": "cells",
                          "total": manifest["admission"]["n_admitted_cells"],
                          "processed": manifest["admission"]["n_admitted_cells"]}],
            "limitations": [*manifest["gaps"], *([] if measurement_path else
                            ["No feature measurement declared: every cell is assumed to have measured every gene."])]}


if __name__ == "__main__":
    main()
