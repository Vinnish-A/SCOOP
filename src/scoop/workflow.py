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

from scsp_agent_sop.config import deep_get, project_root_from_run_root, read_yaml, resolve_run_root

from .io import atomic_json, sha256_file
from .results import archive_receipt, config_sections_sha256, result_payload

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


# run.yaml sections each stage (and the library code it passes the config to)
# reads. A receipt is stale when these change, not when any other section does;
# an undeclared stage binds the whole file.
CONFIG_SECTIONS: dict[str, tuple[str, ...]] = {
    "validate": ("run", "keys", "qc"),
    "qc": ("run", "keys", "qc"),
    "core": ("run", "keys", "qc", "core"),
    "programs": ("run", "keys", "programs"),
    "annotation-evidence": ("run", "annotation"),
    "tumor-cnv": ("run", "keys", "annotation", "tumor_fastcnv"),
    "annotation-review": ("run", "keys", "annotation", "tumor_fastcnv"),
    "annotation-commit": ("run", "keys", "annotation", "tumor_fastcnv"),
    "atlas-integrate": ("run", "keys", "core", "atlas"),
    "atlas-annotate": ("run", "annotation", "atlas"),
    "atlas-report": ("run", "atlas"),
    "atlas-scpoli": (),
    "atlas-validate": (),
    "spatial": ("run", "qc", "spatial"),
    "spatial-build-reference": (),
    "spatial-identifiability": (),
    "ccc": ("run", "ccc"),
    "de-prepare": ("run", "keys", "qc", "de"),
    "abundance-compare": ("abundance",),
    "deliver": ("run", "delivery"),
}


def _default_inputs(stage: str, run_root: Path, config: dict[str, Any]) -> list[Path]:
    """What a command that declares no inputs reads without --input; first existing wins, as in the command."""
    artifacts = run_root / "artifacts"
    source = deep_get(config, "run.input_h5ad")
    source = [project_root_from_run_root(run_root) / source] if source else []
    return {
        "qc": [artifacts / "validated_adata.h5ad", *source],
        "core": [artifacts / "adata_qc.h5ad"],
        "de-prepare": [artifacts / "final_adata.h5ad", artifacts / "adata_annotation_evidence.h5ad"],
        "ccc": [artifacts / "adata_annotation_evidence.h5ad"],
    }.get(stage, [])


def stage_receipt_path(run_root: Path, stage: str, config: dict[str, Any]) -> Path:
    """state/stages/<stage>.json; spatial fits with an output label each keep their own receipt."""
    label = deep_get(config, "spatial.rctd.output_label") if stage == "spatial" else None
    return run_root / "state" / "stages" / (f"{stage}-{label}.json" if label else f"{stage}.json")


def declared_object(run_root: Path, stage: str, config: dict[str, Any]) -> Path | None:
    """The one H5AD or obs sidecar that ``stage``'s receipt last declared, or None without a receipt.

    The receipt decides, not which files happen to exist: an H5AD from an
    earlier run must not displace a newer sidecar.
    """
    receipt_path = stage_receipt_path(run_root, stage, config)
    if not receipt_path.is_file():
        return None
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("status") in {"running", "failed"}:
        raise RuntimeError(f"{stage} is {receipt['status']} in {receipt_path}; rerun it or pass --input.")
    objects = [Path(record["path"]) for record in receipt.get("outputs", [])
               if Path(record["path"]).suffix in {".h5ad", ".parquet"}]
    if len(objects) != 1:
        raise ValueError(f"{receipt_path} declares {len(objects)} H5AD or obs-sidecar outputs, not one; pass --input.")
    return objects[0]


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


def _path_argument(arguments: Iterable[str], flag: str) -> Path | None:
    values = list(arguments)
    for index, value in enumerate(values):
        if value == flag and index + 1 < len(values):
            return Path(values[index + 1]).resolve()
        if value.startswith(f"{flag}="):
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
    config_path = None if help_only else _path_argument(arguments, "--config")
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
        receipt_path = stage_receipt_path(run_root, stage, config)
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
            target = _path_argument(arguments, "--output") or run_root / "artifacts" / artifact_names[stage]
            if not target.is_file():
                status = "failed"
                review_reasons.append(f"Declared stage output is missing: {target}")
            else:
                if not any(Path(record["path"]).resolve() == target.resolve() for record in outputs):
                    outputs.append({"path": str(target), "bytes": target.stat().st_size,
                                    "sha256": sha256_file(target)})
        configuration = {"path": str(config_path), "sha256": config_hash, "role": "configuration"}
        if stage in CONFIG_SECTIONS:
            sections = CONFIG_SECTIONS[stage]
            configuration.update(config_sections=list(sections),
                                 config_sections_sha256=config_sections_sha256(config, sections))
        inputs = [configuration]
        module_path = getattr(module, "__file__", None)
        if module_path and Path(module_path).is_file():
            inputs.append({"path": str(Path(module_path).resolve()), "sha256": sha256_file(module_path), "role": "code"})
        inputs.extend(assessment.get("inputs", []))
        # Commands that return no result still read data; bind it like any declared input.
        data_input = _path_argument(arguments, "--input") or next(
            (path for path in _default_inputs(stage, run_root, config) if path.is_file()), None)
        if (status not in {"skipped", "prepared"} and data_input is not None and data_input.is_file()
                and not any(Path(record["path"]).resolve() == data_input.resolve() for record in inputs)):
            inputs.append({"path": str(data_input.resolve()), "bytes": data_input.stat().st_size,
                           "sha256": sha256_file(data_input), "role": "data"})
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
                "execution": {**assessment.get("execution", {}), "status": execution},
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
