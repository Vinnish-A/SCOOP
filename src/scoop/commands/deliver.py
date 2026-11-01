#!/usr/bin/env python
"""Prune H5AD to minimal delivery state and write manifest.

Usage:
  scoop deliver --config runs/<run_id>/config/run.yaml

This script applies configs/h5ad_schema.yaml, writes final_adata.h5ad and
creates artifacts/manifest.json.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad

from scsp_agent_sop.config import deep_get, project_root_from_run_root, read_yaml, resolve_run_root
from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.manifest import write_manifest
from scsp_agent_sop.storage import ensure_dir, init_file_registry, prune_h5ad, read_adata, register_file, write_h5ad_atomic

from scoop.workflow import declared_object


# The atlas route ends in validated label transfer without per-observation review.
NOT_REVIEWED = ("commit reviewed decisions with `scoop annotate commit` (the atlas route writes no "
                "adata_annotation_committed.h5ad or annotation_reviewed), or export this state explicitly with "
                "--scope current-state")


def review_coverage(adata) -> dict:
    if "annotation_reviewed" not in adata.obs or "annotation_status" not in adata.obs:
        raise ValueError(f"Reviewed annotation delivery requires explicit per-observation review coverage; {NOT_REVIEWED}.")
    flags = adata.obs.annotation_reviewed
    if flags.isna().any() or not flags.isin([True, False, 0, 1]).all() or not flags.all():
        raise ValueError("Unreviewed observations remain; finish their decisions or explicitly export --scope current-state.")
    status = adata.obs.annotation_status
    # excluded_technical is a reviewed outcome (doublet, ambient, ...), not a backlog.
    if not status.isin(["accepted", "labeled", "review_required", "excluded_technical"]).all():
        raise ValueError("Annotation contains unassigned/unreviewed statuses.")
    if "cell_type_lvl3" not in adata.obs:
        raise ValueError("Reviewed annotation delivery requires cell_type_lvl3 labels")
    labels = adata.obs.cell_type_lvl3
    accepted = status.isin(["accepted", "labeled"])
    if (labels[accepted].isna().any() or labels[accepted].astype(str).str.strip().isin(["", "unassigned"]).any()):
        raise ValueError("Accepted observations require assigned labels")
    return {"observations": int(adata.n_obs), "reviewed": int(flags.sum()),
            "accepted": int(status.isin(["accepted", "labeled"]).sum()),
            "reviewed_unresolved": int(status.eq("review_required").sum()),
            "technical_excluded": int(status.eq("excluded_technical").sum())}


def main() -> dict:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--schema", default="configs/h5ad_schema.yaml")
    ap.add_argument("--output", default=None)
    ap.add_argument("--scope", choices=["reviewed-annotation", "current-state"], default=None)
    args = ap.parse_args()
    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    project_root = project_root_from_run_root(run_root)
    scope = args.scope or deep_get(cfg, "delivery.scope", "reviewed-annotation")
    if scope not in {"reviewed-annotation", "current-state"}:
        raise ValueError("delivery.scope must be reviewed-annotation or current-state")
    # annotate commit may have written an obs sidecar; its receipt names what it wrote.
    input_path = (Path(args.input) if args.input else declared_object(run_root, "annotation-commit", cfg)
                  or run_root / "artifacts" / "adata_annotation_committed.h5ad")
    if not input_path.exists() and scope == "current-state" and args.input is None:
        for candidate in [
            "adata_atlas_validated.h5ad",
            "adata_scpoli.h5ad",
            "adata_annotation_committed.h5ad",
            "adata_annotation_evidence.h5ad",
            "adata_programs.h5ad",
            "adata_core.h5ad",
            "adata_qc.h5ad",
        ]:
            c = run_root / "artifacts" / candidate
            if c.exists():
                input_path = c
                break
    if not input_path.is_file():
        raise FileNotFoundError(f"Declared {scope} input is missing: {input_path}."
                                + (f" To deliver, {NOT_REVIEWED}." if scope == "reviewed-annotation" else ""))
    output = (
        Path(args.output)
        if args.output
        else project_root / deep_get(cfg, "run.output_h5ad", str(run_root / "artifacts" / "final_adata.h5ad"))
    )
    schema_path = project_root / args.schema
    if not schema_path.is_file():
        schema_path = Path(__file__).resolve().parents[1] / "data" / "h5ad-schema.yaml"
    schema = read_yaml(schema_path)
    adata = read_adata(input_path)
    coverage = review_coverage(adata) if scope == "reviewed-annotation" else {"observations": int(adata.n_obs)}
    if input_path.resolve() == output.resolve():
        raise ValueError("Delivery output must differ from its input; preserve the source artifact.")
    original_names = adata.obs_names.copy()
    init_file_registry(adata, deep_get(cfg, "run.run_id", run_root.name))
    prune_h5ad(adata, schema)
    if scope == "reviewed-annotation":
        removed = [c for c in ("annotation_reviewed", "annotation_status", "cell_type_lvl3") if c not in adata.obs]
        if removed:
            raise ValueError(f"Delivery schema {schema_path} removed reviewed-annotation fields {removed}; add them to obs_keep.")
        review_coverage(adata)
    # Delivery's own record, like the registry entry below, is not pruned state.
    adata.uns["delivery_scope"] = {"scope": scope, "input": str(input_path.resolve()),
                                    "all_fine_identities_resolved": False, **coverage}
    ensure_dir(output.parent)
    # A file cannot contain its own SHA. The external manifest fingerprints the
    # completed object; the embedded registry contains only its location/schema.
    register_file(adata, key="final_adata", path=output, category="artifacts", schema="h5ad.final_minimal_state.v1",
                  fingerprint=False, extra={"self_reference": True})
    write_h5ad_atomic(adata, output)
    manifest_path = run_root / "artifacts" / "manifest.json"
    log_decision(
        run_root,
        module="prune_and_manifest",
        decision="delivery_h5ad_pruned",
        reason="Final H5AD contains minimal state and registry; large tables remain in run directory.",
        parameters={"schema": args.schema, "output": str(output), "manifest": str(manifest_path)},
    )
    reread = ad.read_h5ad(output, backed="r")
    try:
        if not reread.obs_names.equals(original_names):
            raise ValueError("Delivered observation axis changed")
    finally:
        reread.file.close()
    paths = [output]
    for category in adata.uns["file_registry"].values():
        if isinstance(category, dict):
            paths.extend(record["path"] for record in category.values() if isinstance(record, dict) and "path" in record)
    write_manifest(run_root, paths=paths, metadata={"delivery_scope": scope, "coverage": coverage})
    manifest = json.loads(manifest_path.read_text())
    return {"inputs": [str(input_path.resolve())], "outputs": [str(output.resolve()), str(manifest_path.resolve())],
            "execution": {"status": "completed"},
            "validation": {"status": "limited" if manifest["missing_registered_artifacts"] else "passed",
                           "scope": "delivery_persistence_and_declared_review_coverage",
                           "checks": [{"name": "observation_axis_readback", "status": "passed"}]},
            "interpretation": {"status": "not_assessed"},
            "coverage": [{"scope": scope, "unit": "observation_views", "total": adata.n_obs, "processed": adata.n_obs}],
            "limitations": ["Reviewed unresolved identities remain explicit; export does not establish every fine identity."]}


if __name__ == "__main__":
    main()
