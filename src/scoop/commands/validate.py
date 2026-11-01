#!/usr/bin/env python
"""Validate an input H5AD and initialise SOP state.

Usage:
  scoop validate --config runs/<run_id>/config/run.yaml

This script checks that the input H5AD contains raw counts in the configured
layer, required metadata keys, and unique obs/var names. It writes a validated
H5AD into the run artifacts directory and logs a decision record.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

from scsp_agent_sop.config import deep_get, project_root_from_run_root, read_yaml, resolve_obs_keys, resolve_run_root
from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.storage import init_file_registry, register_file, write_h5ad_atomic


def _require_unique_names(adata) -> None:
    # Renaming here would silently re-key every downstream join to the source.
    for axis, names in (("obs", adata.obs_names), ("var", adata.var_names)):
        duplicated = names[names.duplicated()].unique()
        if len(duplicated):
            raise ValueError(
                f"{axis}_names are not unique: {len(duplicated)} duplicated, e.g. {list(map(str, duplicated[:5]))}. "
                "Resolve the identifiers in the source H5AD."
            )


def _require_integer_counts(matrix, layer: str) -> None:
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).reshape(-1)
    step = 1 << 24  # bounded temporaries on atlas-sized layers
    for start in range(0, values.size, step):
        chunk = values[start : start + step]
        bad = ~(np.isfinite(chunk) & (chunk >= 0) & (chunk == np.floor(chunk)))
        if bad.any():
            raise ValueError(
                f"layers['{layer}'] must hold finite, non-negative integer counts; found {chunk[bad][0]!r}. "
                "Point qc.counts_layer at raw counts."
            )


def main() -> dict:
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
    sample_key = resolve_obs_keys(cfg).sample
    missing = []
    if counts_layer not in adata.layers:
        missing.append(f"layers['{counts_layer}']")
    if sample_key not in adata.obs:
        missing.append(f"obs['{sample_key}']")
    if missing:
        raise ValueError("Missing required fields: " + ", ".join(missing))
    _require_unique_names(adata)
    _require_integer_counts(adata.layers[counts_layer], counts_layer)
    adata.uns["schema_version"] = "scsp_agent_sop.h5ad.v1"
    init_file_registry(adata, run_id=deep_get(cfg, "run.run_id", run_root.name))
    # A file cannot contain its own SHA; the stage receipt fingerprints it.
    register_file(adata, key="validated_adata", path=output, category="artifacts", schema="h5ad.validated.v1",
                  fingerprint=False, extra={"self_reference": True})
    write_h5ad_atomic(adata, output)
    log_decision(
        run_root,
        module="validate_input",
        decision="validated",
        reason="required fields present, obs/var names unique, counts are non-negative integers",
        parameters={"input_h5ad": str(input_path), "output": str(output), "counts_layer": counts_layer},
    )
    return {
        "inputs": [str(input_path.resolve())],
        "outputs": [str(output.resolve())],
        "execution": {"status": "completed"},
        "validation": {
            "status": "passed",
            "scope": "input_h5ad_contract",
            "checks": [
                {"name": "required_fields", "status": "passed", "fields": [f"layers['{counts_layer}']", f"obs['{sample_key}']"]},
                {"name": "unique_obs_var_names", "status": "passed"},
                {"name": "integer_counts", "status": "passed", "layer": counts_layer},
            ],
        },
        "interpretation": {"status": "not_applicable"},
        "coverage": [{"scope": "input_h5ad", "unit": "observations", "total": int(adata.n_obs), "processed": int(adata.n_obs)}],
    }


if __name__ == "__main__":
    main()
