#!/usr/bin/env python
"""Run or prepare RCTD-py spatial deconvolution.

Usage:
  python scripts/05_spatial_rctd.py --config runs/<run_id>/config/run.yaml
  python scripts/05_spatial_rctd.py --config runs/<run_id>/config/run.yaml --dry-run

The default low-resolution mode is `full`. The script builds the configured
RCTD-py command and executes it unless `--dry-run` is set. RCTD output is kept
external; only summaries should be merged into spatial AnnData.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from scsp_agent_sop.config import read_yaml, deep_get, resolve_run_root, project_root_from_run_root
from scsp_agent_sop.spatial import choose_rctd_mode, build_rctd_command, run_rctd_command
from scsp_agent_sop.storage import write_json
from scsp_agent_sop.decision_log import log_decision


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--mode", default=None, choices=["full", "multi", "doublet"])
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    project_root = project_root_from_run_root(run_root)

    if not deep_get(cfg, "spatial.enabled", False):
        print("Spatial module disabled in config. Set spatial.enabled=true to run.")
        return
    spatial_h5ad = project_root / deep_get(cfg, "spatial.spatial_h5ad")
    reference_h5ad = project_root / deep_get(cfg, "spatial.reference_h5ad")
    expected = deep_get(cfg, "spatial.rctd.expected_cells_per_unit")
    mode = args.mode or choose_rctd_mode("spot", expected_cells_per_unit=expected, default_low_res=deep_get(cfg, "spatial.rctd.mode_low_resolution", "full"))
    output_h5ad = run_root / "05_spatial" / "deconvolution" / f"rctd_{mode}.h5ad"
    cmd = build_rctd_command(
        deep_get(cfg, "spatial.rctd.command_template"),
        spatial_h5ad=spatial_h5ad,
        reference_h5ad=reference_h5ad,
        cell_type_key=deep_get(cfg, "spatial.rctd.cell_type_key", "cell_type_lvl3"),
        mode=mode,
        output_h5ad=output_h5ad,
    )
    record = {"mode": mode, "command": cmd, "output_h5ad": str(output_h5ad)}
    write_json(record, run_root / "05_spatial" / "deconvolution" / "rctd_command.json")
    if args.dry_run:
        print(" ".join(cmd))
        return
    completed = run_rctd_command(cmd, cwd=project_root)
    log_decision(
        run_root,
        module="spatial_rctd",
        decision="rctd_command_executed",
        reason="RCTD-py run executed with SOP mode decision.",
        parameters=record,
        evidence={"stdout_tail": completed.stdout[-1000:]},
    )


if __name__ == "__main__":
    main()
