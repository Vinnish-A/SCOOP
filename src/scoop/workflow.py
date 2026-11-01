"""Semantic workflow graph and stage receipt helpers."""

from __future__ import annotations

import importlib
import json
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from scsp_agent_sop.config import read_yaml, resolve_run_root

from .io import atomic_json, sha256_file
from .results import archive_receipt, result_payload

RECEIPT_SCHEMA = "scoop.stage-receipt.v1"

STAGE_DIRECTORIES = {
    "annotation-evidence": "annotation",
    "annotation-review": "annotation",
    "annotation-commit": "annotation",
    "tumor-cnv": "annotation",
    "atlas-integrate": "atlas",
    "atlas-annotate": "atlas",
    "atlas-report": "atlas",
    "atlas-approval": "atlas",
    "atlas-scpoli": "atlas",
    "atlas-validate": "atlas",
    "de-prepare": "differential-expression",
    "de-run": "differential-expression",
    "abundance-compare": "differential-expression",
    "deliver": "delivery",
}


@dataclass(frozen=True)
class StageSpec:
    name: str
    depends_on: tuple[str, ...] = ()
    capability: str = "core"
    review_gate: bool = False


WORKFLOWS: dict[str, tuple[StageSpec, ...]] = {
    "abundance": (StageSpec("abundance-compare"),),
    "standard": (
        StageSpec("validate"),
        StageSpec("qc", ("validate",)),
        StageSpec("core", ("qc",)),
        StageSpec("programs", ("core",)),
        StageSpec("annotation-evidence", ("programs",)),
        StageSpec("annotation-review", ("annotation-evidence",), review_gate=True),
        StageSpec("annotation-commit", ("annotation-review",)),
        StageSpec("deliver", ("annotation-commit",)),
    ),
    "tumor": (
        StageSpec("validate"),
        StageSpec("qc", ("validate",)),
        StageSpec("core", ("qc",)),
        StageSpec("programs", ("core",)),
        StageSpec("annotation-evidence", ("programs",)),
        StageSpec("tumor-cnv", ("annotation-evidence",)),
        StageSpec("annotation-review", ("tumor-cnv",), review_gate=True),
        StageSpec("annotation-commit", ("annotation-review",)),
        StageSpec("deliver", ("annotation-commit",)),
    ),
    "atlas": (
        StageSpec("validate"),
        StageSpec("qc", ("validate",)),
        StageSpec("atlas-integrate", ("qc",), capability="atlas"),
        StageSpec("programs", ("atlas-integrate",)),
        StageSpec("atlas-annotate", ("programs",), capability="atlas"),
        StageSpec("atlas-report", ("atlas-annotate",), capability="atlas", review_gate=True),
        StageSpec("atlas-approval", ("atlas-report",), capability="atlas", review_gate=True),
        StageSpec("atlas-scpoli", ("atlas-approval",), capability="atlas"),
        StageSpec("atlas-validate", ("atlas-scpoli",), capability="atlas"),
        StageSpec("deliver", ("atlas-validate",)),
    ),
    "spatial": (
        StageSpec("annotation-commit"),
        StageSpec("spatial", ("annotation-commit",), capability="spatial"),
    ),
    "ccc": (
        StageSpec("annotation-commit"),
        StageSpec("ccc", ("annotation-commit",), capability="ccc"),
    ),
    "differential-expression": (
        StageSpec("annotation-commit"),
        StageSpec("de-prepare", ("annotation-commit",)),
        StageSpec("de-run", ("de-prepare",)),
    ),
}


def workflow_plan(name: str, config_path: str | Path) -> dict[str, Any]:
    if name not in WORKFLOWS:
        raise ValueError(f"unknown workflow: {name}")
    config_path = Path(config_path).resolve()
    config = read_yaml(config_path)
    run_root = resolve_run_root(config_path, config).resolve()
    return {
        "schema_version": "scoop.workflow-plan.v2",
        "workflow": name,
        "run_id": str(config.get("run", {}).get("run_id", run_root.name)),
        "config": str(config_path),
        "run_root": str(run_root),
        "stages": [asdict(stage) for stage in WORKFLOWS[name]],
    }


def _config_argument(arguments: Iterable[str]) -> Path | None:
    values = list(arguments)
    for index, value in enumerate(values):
        if value == "--config" and index + 1 < len(values):
            return Path(values[index + 1]).resolve()
        if value.startswith("--config="):
            return Path(value.split("=", 1)[1]).resolve()
    return None


def _artifact_records(root: Path) -> list[dict[str, Any]]:
    if not root.exists():
        return []
    records: list[dict[str, Any]] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        record: dict[str, Any] = {
            "path": str(path),
            "bytes": path.stat().st_size,
        }
        if path.stat().st_size <= 64 * 1024 * 1024:
            record["sha256"] = sha256_file(path)
        else:
            record["sha256"] = None
            record["hash_status"] = "deferred-large-artifact"
        records.append(record)
    return records


def run_command_module(stage: str, module_name: str, arguments: list[str]) -> int:
    help_only = any(value in {"-h", "--help"} for value in arguments)
    config_path = None if help_only else _config_argument(arguments)
    run_root: Path | None = None
    receipt_path: Path | None = None
    started = time.monotonic()
    started_at = datetime.now(timezone.utc).isoformat()
    config_hash: str | None = None
    decision_offset = 0
    config: dict[str, Any] = {}
    if config_path and config_path.is_file():
        config = read_yaml(config_path)
        run_root = resolve_run_root(config_path, config).resolve()
        receipt_path = run_root / "state" / "stages" / f"{stage}.json"
        config_hash = sha256_file(config_path)
        decision_log = run_root / "logs/decision_log.jsonl"
        decision_offset = decision_log.stat().st_size if decision_log.is_file() else 0
        archive_receipt(receipt_path)
        atomic_json(
            receipt_path,
            {
                "schema_version": RECEIPT_SCHEMA,
                "stage": stage,
                "status": "running",
                "config": str(config_path),
                "config_hash": config_hash,
                "started_at": started_at,
            },
        )

    previous = sys.argv
    sys.argv = [f"scoop {stage}", *arguments]
    try:
        module = importlib.import_module(module_name)
        try:
            result = module.main()
        except SystemExit as exc:
            if help_only and exc.code in {None, 0}:
                return 0
            raise
        assessment = result_payload(result, base=config_path.parent if config_path else Path.cwd(), route="managed") if isinstance(result, dict) else {}
        return_code = int(assessment.get("exit_code", 0) if assessment else (result or 0))
        if return_code:
            raise RuntimeError(f"stage returned non-zero status: {return_code}")
    except BaseException as exc:
        if receipt_path is not None:
            atomic_json(
                receipt_path,
                {
                    "schema_version": RECEIPT_SCHEMA,
                    "stage": stage,
                    "status": "failed",
                    "config": str(config_path),
                    "config_hash": config_hash,
                    "error": {"type": type(exc).__name__, "message": str(exc)},
                    "started_at": started_at,
                    "finished_at": datetime.now(timezone.utc).isoformat(),
                    "elapsed_seconds": round(time.monotonic() - started, 6),
                },
            )
        raise
    finally:
        sys.argv = previous

    if receipt_path is not None and run_root is not None:
        review_required = stage in {"annotation-review", "atlas-report", "atlas-approval"}
        review_reasons = []
        decision_log = run_root / "logs/decision_log.jsonl"
        if decision_log.is_file():
            with decision_log.open("rb") as handle:
                handle.seek(decision_offset)
                for line in handle:
                    entry = json.loads(line)
                    if entry.get("human_review_required"):
                        review_required = True
                        review_reasons.append(entry.get("review_reason") or entry["reason"])
        status = "awaiting_review" if review_required else "executed"
        if stage in {"spatial", "ccc"} and not config.get(stage, {}).get("enabled", False):
            status = "skipped"
        elif "--dry-run" in arguments:
            status = "prepared"
        # A directory inventory can contain old files. Only declared outputs
        # belong to this execution; persistence alone is not numerical validation.
        outputs = list(assessment.get("outputs", []))
        artifact_names = {"qc": "adata_qc.h5ad", "core": "adata_core.h5ad", "programs": "adata_programs.h5ad"}
        if stage in artifact_names and status in {"executed", "awaiting_review"}:
            target = run_root / "artifacts" / artifact_names[stage]
            for index, value in enumerate(arguments):
                if value == "--output" and index + 1 < len(arguments):
                    target = Path(arguments[index + 1]).resolve()
                elif value.startswith("--output="):
                    target = Path(value.split("=", 1)[1]).resolve()
            if not target.is_file():
                status = "failed"
                review_reasons.append(f"Declared stage output is missing: {target}")
            else:
                if not any(Path(record["path"]).resolve() == target.resolve() for record in outputs):
                    outputs.append({"path": str(target), "bytes": target.stat().st_size,
                                    "sha256": sha256_file(target)})
        inputs = [{"path": str(config_path), "sha256": config_hash, "role": "configuration"}]
        module_path = getattr(module, "__file__", None)
        if module_path and Path(module_path).is_file():
            inputs.append({"path": str(Path(module_path).resolve()), "sha256": sha256_file(module_path), "role": "code"})
        inputs.extend(assessment.get("inputs", []))
        execution = "not_run" if status in {"prepared", "skipped"} else "failed" if status == "failed" else "completed"
        if assessment.get("execution", {}).get("status") == "failed":
            status, execution = "failed", "failed"
        elif assessment.get("execution", {}).get("status") == "not_run" and status != "failed":
            status, execution = "skipped", "not_run"
        atomic_json(
            receipt_path,
            {
                **assessment,
                "schema_version": RECEIPT_SCHEMA,
                "stage": stage,
                "status": status,
                "config": str(config_path),
                "config_hash": config_hash,
                "outputs": outputs,
                "inputs": inputs,
                "execution": {"status": execution},
                "validation": assessment.get("validation", {"status": "not_assessed"}),
                "interpretation": assessment.get("interpretation", {"status": "not_assessed"}),
                "coverage": assessment.get("coverage", []),
                "route": "managed",
                "review_required": review_required,
                "review_reason": "; ".join(review_reasons) if review_reasons else ("decision required" if review_required else None),
                "started_at": started_at,
                "finished_at": datetime.now(timezone.utc).isoformat(),
                "elapsed_seconds": round(time.monotonic() - started, 6),
            },
        )
    return 1 if receipt_path is not None and run_root is not None and status == "failed" else 0
