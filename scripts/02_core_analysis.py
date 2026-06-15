#!/usr/bin/env python
"""Run the 02_core analysis through FastCore.

Usage:
  python scripts/02_core_analysis.py --config runs/<run_id>/config/run.yaml

The SOP entry point remains 02_core. FastCore selects one backend before
execution and uses Scanpy legacy only as the single whole-pipeline fallback.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import anndata as ad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fastcore.backend_plan import plan_fastcore_backend
from scsp_agent_sop.config import read_yaml, deep_get, resolve_run_root
from scsp_agent_sop.core_runner import run_core_pipeline
from scsp_agent_sop.storage import init_file_registry, ensure_dir


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--output", default=None)
    args = ap.parse_args()
    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    input_path = Path(args.input) if args.input else run_root / "artifacts" / "adata_qc.h5ad"
    output = Path(args.output) if args.output else run_root / "artifacts" / "adata_core.h5ad"
    preplan = plan_fastcore_backend(cfg, input_path=input_path) if deep_get(cfg, "core.engine", "fastcore") == "fastcore" else None
    if preplan is not None and preplan.selected_backend == "omicverse_rust_oom":
        adata = None
    else:
        adata = ad.read_h5ad(input_path)
        init_file_registry(adata, deep_get(cfg, "run.run_id", run_root.name))

    run_core_pipeline(
        adata,
        cfg,
        run_root,
        input_path=input_path,
        output_path=output,
    )

    if adata is not None:
        ensure_dir(output.parent)
        adata.write_h5ad(output)


if __name__ == "__main__":
    main()
