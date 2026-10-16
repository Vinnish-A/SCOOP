"""Unified semantic CLI for SCOOP."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Sequence

import yaml

from scsp_agent_sop import __version__

from . import approval
from .doctor import emit as emit_doctor
from .io import atomic_json
from .release import rollback as rollback_release
from .release import status as release_status
from .workflow import WORKFLOWS, run_command_module, workflow_plan

COMMANDS: dict[tuple[str, ...], tuple[str, str]] = {
    ("validate",): ("validate", "scoop.commands.validate"),
    ("qc",): ("qc", "scoop.commands.qc"),
    ("core",): ("core", "scoop.commands.core"),
    ("programs",): ("programs", "scoop.commands.programs"),
    ("atlas", "integrate"): ("atlas-integrate", "scoop.commands.atlas_integrate"),
    ("atlas", "annotate"): ("atlas-annotate", "scoop.commands.atlas_annotate"),
    ("atlas", "report"): ("atlas-report", "scoop.commands.atlas_report"),
    ("atlas", "scpoli"): ("atlas-scpoli", "scoop.commands.atlas_scpoli"),
    ("atlas", "validate"): ("atlas-validate", "scoop.commands.atlas_validate"),
    ("annotate", "evidence"): ("annotation-evidence", "scoop.commands.annotation_evidence"),
    ("annotate", "tumor-cnv"): ("tumor-cnv", "scoop.commands.tumor_cnv"),
    ("annotate", "review"): ("annotation-review", "scoop.commands.annotation_review"),
    ("annotate", "commit"): ("annotation-commit", "scoop.commands.annotation_commit"),
    ("spatial",): ("spatial", "scoop.commands.spatial"),
    ("ccc",): ("ccc", "scoop.commands.ccc"),
    ("de", "prepare"): ("de-prepare", "scoop.commands.de_prepare"),
    ("de", "run"): ("de-run", "scoop.commands.de_run"),
    ("abundance", "compare"): ("abundance-compare", "scoop.commands.abundance_compare"),
    ("deliver",): ("deliver", "scoop.commands.deliver"),
}


def _route(argv: list[str]) -> tuple[str, str, list[str]] | None:
    for key in sorted(COMMANDS, key=len, reverse=True):
        if tuple(argv[: len(key)]) == key:
            stage, module = COMMANDS[key]
            return stage, module, argv[len(key) :]
    return None


def _data_path(name: str) -> Path:
    return Path(__file__).resolve().parent / "data" / name


def _cmd_init(args: argparse.Namespace) -> int:
    project = Path(args.project).expanduser().resolve()
    run_root = project / "runs" / args.run_id
    config_dir = run_root / "config"
    config_dir.mkdir(parents=True, exist_ok=True)
    for name in (
        "validate",
        "qc",
        "core",
        "programs",
        "annotation",
        "atlas",
        "spatial",
        "ccc",
        "differential-expression",
        "delivery",
        "artifacts",
        "reports",
        "logs",
        "state/stages",
    ):
        (run_root / name).mkdir(parents=True, exist_ok=True)
    config = yaml.safe_load(_data_path("default-run.yaml").read_text(encoding="utf-8"))
    config.setdefault("run", {})["run_id"] = args.run_id
    config["run"]["input_h5ad"] = args.input_h5ad
    config["run"]["output_h5ad"] = str(run_root / "artifacts" / "final_adata.h5ad")
    target = config_dir / "run.yaml"
    target.write_text(yaml.safe_dump(config, sort_keys=False, allow_unicode=True), encoding="utf-8")
    print(target)
    return 0


def _cmd_plan(args: argparse.Namespace) -> int:
    payload = workflow_plan(args.workflow, args.config)
    run_root = Path(payload["run_root"])
    output = Path(args.output).resolve() if args.output else run_root / "state" / "workflow-plan.json"
    atomic_json(output, payload)
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


def _cmd_atlas(args: argparse.Namespace) -> int:
    if args.atlas_action == "approve":
        path = approval.approve(
            gate_path=args.gate,
            run_id=args.run_id,
            annotation_version=args.annotation_version,
            approved_by=args.approved_by,
            input_h5ad=args.input_h5ad,
            annotation_report=args.annotation_report,
            label_table=args.label_table,
            labeled_indices=args.labeled_indices,
            notes=args.notes,
        )
        gate = Path(path).resolve()
        run_root = gate.parent.parent if gate.parent.name == "atlas" else None
        if run_root is not None:
            verified = approval.verify(gate)
            atomic_json(
                run_root / "state" / "stages" / "atlas-approval.json",
                {
                    "schema_version": "scoop.stage-receipt.v1",
                    "stage": "atlas-approval",
                    "status": "committed" if verified["ok"] else "failed",
                    "review_required": False,
                    "approval_gate": str(gate),
                    "annotation_version": args.annotation_version,
                    "approved_by": args.approved_by,
                    "errors": verified["errors"],
                },
            )
        print(path)
        return 0
    result = approval.verify(args.gate)
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if result["ok"] else 2


def _cmd_runtime(args: argparse.Namespace) -> int:
    if args.runtime_action == "rollback":
        print(rollback_release())
        return 0
    payload = release_status()
    print(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if payload["ok"] else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="scoop", description="Single Cell Omics Operating Protocol")
    parser.add_argument("--version", action="version", version=f"SCOOP {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    doctor = commands.add_parser("doctor", help="inspect managed runtime capabilities")
    doctor.add_argument("--json", action="store_true")
    doctor.set_defaults(func=lambda args: emit_doctor(json_output=args.json))

    init = commands.add_parser("init", help="create a semantic project run")
    init.add_argument("--project", required=True)
    init.add_argument("--run-id", required=True)
    init.add_argument("--input-h5ad", required=True)
    init.set_defaults(func=_cmd_init)

    plan = commands.add_parser("plan", help="write a semantic workflow plan")
    plan.add_argument("--workflow", required=True, choices=sorted(WORKFLOWS))
    plan.add_argument("--config", required=True)
    plan.add_argument("--output")
    plan.set_defaults(func=_cmd_plan)

    atlas = commands.add_parser("atlas", help="manage the scVI annotation approval gate")
    atlas_actions = atlas.add_subparsers(dest="atlas_action", required=True)
    approve_parser = atlas_actions.add_parser("approve")
    approve_parser.add_argument("--gate", required=True)
    approve_parser.add_argument("--run-id", required=True)
    approve_parser.add_argument("--annotation-version", required=True)
    approve_parser.add_argument("--approved-by", required=True)
    approve_parser.add_argument("--input-h5ad", required=True)
    approve_parser.add_argument("--annotation-report", required=True)
    approve_parser.add_argument("--label-table", required=True)
    approve_parser.add_argument("--labeled-indices", required=True)
    approve_parser.add_argument("--notes")
    approve_parser.set_defaults(func=_cmd_atlas)
    verify_parser = atlas_actions.add_parser("verify-gate")
    verify_parser.add_argument("--gate", required=True)
    verify_parser.set_defaults(func=_cmd_atlas)

    runtime = commands.add_parser("runtime", help="inspect or roll back the managed release")
    runtime_actions = runtime.add_subparsers(dest="runtime_action", required=True)
    runtime_status = runtime_actions.add_parser("status")
    runtime_status.set_defaults(func=_cmd_runtime)
    runtime_rollback = runtime_actions.add_parser("rollback")
    runtime_rollback.set_defaults(func=_cmd_runtime)
    return parser


def _print_stage_help() -> None:
    print("Semantic analysis commands:")
    for key in sorted(COMMANDS):
        print("  scoop " + " ".join(key) + " ...")


def main(argv: Sequence[str] | None = None) -> int:
    values = list(argv if argv is not None else sys.argv[1:])
    routed = _route(values)
    if routed is not None:
        stage, module, remaining = routed
        return run_command_module(stage, module, remaining)
    if not values:
        parser = build_parser()
        parser.print_help()
        _print_stage_help()
        return 0
    try:
        args = build_parser().parse_args(values)
        return int(args.func(args))
    except KeyboardInterrupt:
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
