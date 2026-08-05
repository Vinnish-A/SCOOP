#!/usr/bin/env python
"""Prune H5AD to minimal delivery state and write manifest.

Usage:
  scoop deliver --config runs/<run_id>/config/run.yaml

This script applies configs/h5ad_schema.yaml, writes final_adata.h5ad and
creates artifacts/manifest.json.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad

from scsp_agent_sop.config import deep_get, project_root_from_run_root, read_yaml, resolve_run_root
from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.manifest import write_manifest
from scsp_agent_sop.storage import ensure_dir, init_file_registry, prune_h5ad, register_file


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--schema", default="configs/h5ad_schema.yaml")
    ap.add_argument("--output", default=None)
    args = ap.parse_args()
    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    project_root = project_root_from_run_root(run_root)
    input_path = Path(args.input) if args.input else run_root / "artifacts" / "adata_ccc_indexed.h5ad"
    if not input_path.exists():
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
    output = (
        Path(args.output)
        if args.output
        else project_root / deep_get(cfg, "run.output_h5ad", str(run_root / "artifacts" / "final_adata.h5ad"))
    )
    schema_path = project_root / args.schema
    if not schema_path.is_file():
        schema_path = Path(__file__).resolve().parents[1] / "data" / "h5ad-schema.yaml"
    schema = read_yaml(schema_path)
    adata = ad.read_h5ad(input_path)
    init_file_registry(adata, deep_get(cfg, "run.run_id", run_root.name))
    prune_h5ad(adata, schema)
    ensure_dir(output.parent)
    adata.write_h5ad(output)
    register_file(adata, key="final_adata", path=output, category="artifacts", schema="h5ad.final_minimal_state.v1")
    adata.write_h5ad(output)
    manifest_path = run_root / "artifacts" / "manifest.json"
    log_decision(
        run_root,
        module="prune_and_manifest",
        decision="delivery_h5ad_pruned",
        reason="Final H5AD contains minimal state and registry; large tables remain in run directory.",
        parameters={"schema": args.schema, "output": str(output), "manifest": str(manifest_path)},
    )
    write_manifest(run_root)


if __name__ == "__main__":
    main()
