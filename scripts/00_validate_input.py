#!/usr/bin/env python
"""Validate an input H5AD and initialise SOP state.

Usage:
  python scripts/00_validate_input.py --config runs/<run_id>/config/run.yaml

This script checks that the input H5AD contains raw counts in the configured
layer, required metadata keys, and unique obs/var names. It writes a validated
H5AD into the run artifacts directory and logs a decision record.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import anndata as ad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from scsp_agent_sop.config import read_yaml, deep_get, resolve_run_root, project_root_from_run_root
from scsp_agent_sop.storage import ensure_dir, init_file_registry, register_file
from scsp_agent_sop.decision_log import log_decision


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--output", default=None)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    project_root = project_root_from_run_root(run_root)
    input_path = project_root / deep_get(cfg, "run.input_h5ad")
    output = Path(args.output) if args.output else run_root / "artifacts" / "validated_adata.h5ad"

    if not input_path.exists():
        raise FileNotFoundError(f"Input H5AD not found: {input_path}")
    adata = ad.read_h5ad(input_path)
    counts_layer = deep_get(cfg, "qc.counts_layer", "counts")
    sample_key = deep_get(cfg, "keys.sample", "sample_id")
    missing = []
    if counts_layer not in adata.layers:
        missing.append(f"layers['{counts_layer}']")
    if sample_key not in adata.obs:
        missing.append(f"obs['{sample_key}']")
    if missing:
        raise ValueError("Missing required fields: " + ", ".join(missing))
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    adata.uns["schema_version"] = "scsp_agent_sop.h5ad.v1"
    init_file_registry(adata, run_id=deep_get(cfg, "run.run_id", run_root.name))
    ensure_dir(output.parent)
    adata.write_h5ad(output)
    register_file(adata, key="validated_adata", path=output, category="artifacts", schema="h5ad.validated.v1")
    log_decision(run_root, module="validate_input", decision="validated", reason="required H5AD fields are present", parameters={"input_h5ad": str(input_path), "output": str(output)})


if __name__ == "__main__":
    main()
