"""Semantic workflow graph and stage receipt helpers."""

from __future__ import annotations

import importlib
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from scsp_agent_sop.config import read_yaml, resolve_run_root

from .io import atomic_json, sha256_file

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
    "deliver": "delivery",
}


@dataclass(frozen=True)
class StageSpec:
    name: str
    depends_on: tuple[str, ...] = ()
    profile: str = "core"
    review_gate: bool = False


WORKFLOWS: dict[str, tuple[StageSpec, ...]] = {
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
        StageSpec("atlas-integrate", ("qc",), profile="atlas"),
        StageSpec("programs", ("atlas-integrate",)),
        StageSpec("atlas-annotate", ("programs",), profile="atlas"),
        StageSpec("atlas-report", ("atlas-annotate",), profile="atlas", review_gate=True),
        StageSpec("atlas-approval", ("atlas-report",), profile="atlas", review_gate=True),
        StageSpec("atlas-scpoli", ("atlas-approval",), profile="atlas"),
        StageSpec("atlas-validate", ("atlas-scpoli",), profile="atlas"),
        StageSpec("deliver", ("atlas-validate",)),
    ),
    "spatial": (
        StageSpec("annotation-commit"),
        StageSpec("spatial", ("annotation-commit",), profile="spatial"),
    ),
    "ccc": (
        StageSpec("annotation-commit"),
        StageSpec("ccc", ("annotation-commit",), profile="ccc"),
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
        "schema_version": "scoop.workflow-plan.v1",
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
    config_path = _config_argument(arguments)
    run_root: Path | None = None
    receipt_path: Path | None = None
    started = time.monotonic()
    started_at = datetime.now(timezone.utc).isoformat()
    config_hash: str | None = None
    if config_path and config_path.is_file():
        config = read_yaml(config_path)
        run_root = resolve_run_root(config_path, config).resolve()
        receipt_path = run_root / "state" / "stages" / f"{stage}.json"
        config_hash = sha256_file(config_path)
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
        result = module.main()
        return_code = int(result or 0)
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
        stage_root = run_root / STAGE_DIRECTORIES.get(stage, stage.split("-", 1)[0])
        review_required = stage in {"annotation-review", "atlas-report", "atlas-approval"}
        atomic_json(
            receipt_path,
            {
                "schema_version": RECEIPT_SCHEMA,
                "stage": stage,
                "status": "awaiting_review" if review_required else "committed",
                "config": str(config_path),
                "config_hash": config_hash,
                "outputs": _artifact_records(stage_root),
                "review_required": review_required,
                "review_reason": "human decision required" if review_required else None,
                "started_at": started_at,
                "finished_at": datetime.now(timezone.utc).isoformat(),
                "elapsed_seconds": round(time.monotonic() - started, 6),
            },
        )
    return 0
