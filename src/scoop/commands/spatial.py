#!/usr/bin/env python
"""Run or prepare RCTD-py spatial deconvolution.

Usage:
  scoop spatial --config runs/<run_id>/config/run.yaml
  scoop spatial --config runs/<run_id>/config/run.yaml --dry-run

The default low-resolution mode is `full`. The command builds the configured
RCTD-py command and executes it unless `--dry-run` is set. RCTD output is kept
external; only summaries should be merged into spatial AnnData.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from scsp_agent_sop.config import deep_get, project_root_from_run_root, read_yaml, resolve_run_root
from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.rctd_adapter import run_rctd_api
from scsp_agent_sop.spatial import build_rctd_command, choose_rctd_mode, run_rctd_command
from scsp_agent_sop.storage import write_json


def main() -> dict | None:
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
    reference = deep_get(cfg, "spatial.reference_h5ad")
    reference_h5ad = project_root / reference if reference else None
    options = cfg["spatial"].setdefault("rctd", {})
    for key in ("reference_profiles", "reference_manifest", "spatial_registry", "reference_registry", "q_matrices"):
        if options.get(key):
            options[key] = str((project_root / Path(options[key]).expanduser()).resolve())
    if reference_h5ad is None and not options.get("reference_profiles"):
        raise ValueError("Provide spatial.reference_h5ad or reviewed reference_profiles + reference_manifest")
    expected = deep_get(cfg, "spatial.rctd.expected_cells_per_unit")
    mode = args.mode or choose_rctd_mode(
        deep_get(cfg, "spatial.unit_type", "spot"),
        expected_cells_per_unit=expected,
        default_low_res=deep_get(cfg, "spatial.rctd.mode_low_resolution", "full"),
    )
    output_h5ad = run_root / "spatial" / "deconvolution" / f"rctd_{mode}.h5ad"
    adapter = deep_get(cfg, "spatial.rctd.adapter", "api")
    if adapter not in {"api", "external"}:
        raise ValueError("spatial.rctd.adapter must be api or external.")
    record = {
        "mode": mode,
        "adapter": adapter,
        "output_h5ad": str(output_h5ad),
        "spatial_h5ad": str(spatial_h5ad),
        "reference_h5ad": str(reference_h5ad) if reference_h5ad else None,
        "reference_profiles": options.get("reference_profiles"),
    }
    if adapter == "api":
        write_json(record, run_root / "spatial/deconvolution/rctd_command.json")
        if args.dry_run:
            print(record)
            return
        info = run_rctd_api(spatial_h5ad, reference_h5ad, output_h5ad, cfg, mode)
        log_decision(
            run_root,
            module="spatial_rctd",
            decision="rctd_api_complete",
            reason="Pinned RCTD API executed; fitted and visualization weights kept separately.",
            parameters=record,
            evidence=info,
        )
        return {"inputs": info["inputs"], "outputs": info["outputs"],
                "execution": {"status": "completed"},
                "validation": {"status": "limited", "scope": "native_counts_source_masks_fit_and_persistence",
                               "checks": [{"name": "finite_fitted_weights_exact_component_axis", "status": "passed"},
                                          {"name": "saved_original_observation_axis", "status": "passed"}]},
                "interpretation": {"status": "not_assessed"},
                "coverage": [{"scope": "spatial_fit_source_observations", "unit": "observation_views",
                              "total": info["observations_input"], "processed": info["observations_fitted"]}],
                "limitations": ["RNA allocation weights require reference coverage, identifiability and biological validation.",
                                *info["reference"]["gaps"]],
                "computation_checkpoint": info["computation_checkpoint"]}
    if reference_h5ad is None:
        raise ValueError("The external command adapter requires reference_h5ad")
    cmd = build_rctd_command(
        deep_get(cfg, "spatial.rctd.command_template"),
        spatial_h5ad=spatial_h5ad,
        reference_h5ad=reference_h5ad,
        cell_type_key=deep_get(cfg, "spatial.rctd.cell_type_key", "cell_type_lvl3"),
        mode=mode,
        output_h5ad=output_h5ad,
    )
    record["command"] = cmd
    write_json(record, run_root / "spatial/deconvolution/rctd_command.json")
    if args.dry_run:
        print(" ".join(cmd))
        return
    completed = run_rctd_command(cmd, cwd=project_root)
    if not output_h5ad.is_file():
        raise RuntimeError("External RCTD exited successfully without its declared H5AD output.")
    log_decision(
        run_root,
        module="spatial_rctd",
        decision="rctd_command_executed",
        reason="RCTD-py run executed with SOP mode decision.",
        parameters=record,
        evidence={"stdout_tail": completed.stdout[-1000:]},
    )
    return {"inputs": [str(spatial_h5ad.resolve()), str(reference_h5ad.resolve())],
            "outputs": [str(output_h5ad.resolve())], "execution": {"status": "completed"},
            "validation": {"status": "not_assessed"}, "interpretation": {"status": "not_assessed"},
            "limitations": ["External command returned successfully; its numerical output has not been assessed here."]}


if __name__ == "__main__":
    main()
