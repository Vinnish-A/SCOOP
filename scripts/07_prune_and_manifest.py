#!/usr/bin/env python
"""Prune H5AD to minimal delivery state and write manifest.

Usage:
  python scripts/07_prune_and_manifest.py --config runs/<run_id>/config/run.yaml

This script applies configs/h5ad_schema.yaml, writes final_adata.h5ad and
creates artifacts/manifest.json.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import anndata as ad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from scsp_agent_sop.config import read_yaml, deep_get, resolve_run_root, project_root_from_run_root
from scsp_agent_sop.storage import prune_h5ad, ensure_dir, register_file, init_file_registry
from scsp_agent_sop.manifest import write_manifest
from scsp_agent_sop.decision_log import log_decision


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
        for candidate in ["adata_annotation_evidence.h5ad", "adata_programs.h5ad", "adata_core.h5ad", "adata_qc.h5ad"]:
            c = run_root / "artifacts" / candidate
            if c.exists():
                input_path = c
                break
    output = Path(args.output) if args.output else project_root / deep_get(cfg, "run.output_h5ad", str(run_root / "artifacts" / "final_adata.h5ad"))
    schema = read_yaml(project_root / args.schema)
    adata = ad.read_h5ad(input_path)
    init_file_registry(adata, deep_get(cfg, "run.run_id", run_root.name))
    prune_h5ad(adata, schema)
    ensure_dir(output.parent)
    adata.write_h5ad(output)
    register_file(adata, key="final_adata", path=output, category="artifacts", schema="h5ad.final_minimal_state.v1")
    manifest_path = write_manifest(run_root)
    log_decision(run_root, module="prune_and_manifest", decision="delivery_h5ad_pruned", reason="Final H5AD contains minimal state and registry; large tables remain in run directory.", parameters={"schema": args.schema, "output": str(output), "manifest": str(manifest_path)})


if __name__ == "__main__":
    main()
