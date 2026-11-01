"""Hash-bound results with separate execution, validation and interpretation."""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

from .io import atomic_json, sha256_file

RESULT_SCHEMA = "scoop.result.v1"
ASSESSMENT_STATES = {"passed", "failed", "limited", "not_assessed", "not_applicable"}


def archive_receipt(path: Path) -> None:
    if path.is_file():
        target = path.parent / "history" / path.stem / f"{time.time_ns()}.json"
        atomic_json(target, json.loads(path.read_text(encoding="utf-8")))


def _files(values: list, base: Path) -> list[dict[str, Any]]:
    records = []
    for value in values:
        item = {"path": value} if isinstance(value, str) else dict(value)
        path = Path(item["path"]).expanduser()
        path = (base / path).resolve() if not path.is_absolute() else path.resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Declared result artifact is missing: {path}")
        actual = sha256_file(path)
        if item.get("sha256") and actual != item["sha256"]:
            raise ValueError(f"Result artifact hash differs: {path}")
        records.append({**item, "path": str(path), "bytes": path.stat().st_size, "sha256": actual})
    return records


def result_payload(payload: dict, *, base: Path, route: str) -> dict[str, Any]:
    """Bind declarations to files; artifact integrity does not adjudicate biology."""
    result = dict(payload)
    execution = result.get("execution", {"status": "completed"})
    if execution.get("status") not in {"completed", "failed", "not_run"}:
        raise ValueError("execution.status must be completed, failed or not_run")
    validation = result.get("validation", {"status": "not_assessed", "checks": []})
    if validation.get("status") not in ASSESSMENT_STATES:
        raise ValueError("Unsupported validation.status")
    checks = validation.get("checks", [])
    if validation["status"] == "passed" and (not checks or any(c.get("status") != "passed" for c in checks)):
        raise ValueError("Passed validation requires nonempty, passed checks")
    interpretation = result.get("interpretation", {"status": "not_assessed"})
    if interpretation.get("status") not in {"supported", "unsupported", "not_assessed", "not_applicable"}:
        raise ValueError("Unsupported interpretation.status")
    if interpretation["status"] == "supported" and not interpretation.get("scope"):
        raise ValueError("Supported interpretation requires an explicit scope")
    coverage = result.get("coverage", [])
    for row in coverage:
        if not row.get("unit") or not row.get("scope"):
            raise ValueError("Coverage needs a measured unit and named scope")
        total, processed = row["total"], row["processed"]
        if type(total) is not int or type(processed) is not int or not 0 <= processed <= total:
            raise ValueError("Coverage must satisfy 0 <= processed <= total")
    inputs = _files(result.get("inputs", []), base)
    outputs = _files(result.get("outputs", []), base)
    if execution["status"] == "completed" and not outputs:
        raise ValueError("A registered completed result requires actual output artifacts")
    if validation["status"] == "passed" and not inputs:
        raise ValueError("Passed validation requires bound inputs")
    return {**result, "schema_version": RESULT_SCHEMA, "route": route,
            "execution": execution, "validation": validation, "interpretation": interpretation,
            "coverage": coverage, "inputs": inputs, "outputs": outputs}


def register_result(run_root: Path, manifest: Path) -> Path:
    manifest = manifest.resolve()
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    stage = payload["stage"]
    if not isinstance(stage, str) or not stage or any(c not in "abcdefghijklmnopqrstuvwxyz0123456789-" for c in stage):
        raise ValueError("stage must be a lowercase stage identifier")
    result = result_payload(payload, base=manifest.parent, route="project_adapter")
    result["source_manifest"] = {"path": str(manifest), "sha256": sha256_file(manifest)}
    result["status"] = "failed" if result["execution"]["status"] == "failed" else "executed"
    target = run_root / "state/stages" / f"{stage}.json"
    archive_receipt(target)
    return atomic_json(target, result)


def _assessment(receipt: dict) -> dict:
    if "execution" in receipt:
        defaults = {"validation": {"status": "not_assessed"}, "interpretation": {"status": "not_assessed"},
                    "coverage": [], "limitations": [], "next_actions": []}
        return {"execution": receipt["execution"], **{k: receipt.get(k) or v for k, v in defaults.items()}}
    old = receipt.get("status")
    execution = "completed" if old in {"committed", "awaiting_review", "executed"} else "failed" if old == "failed" else "not_run"
    return {"execution": {"status": execution}, "validation": {"status": "not_assessed"},
            "interpretation": {"status": "not_assessed"}, "coverage": [],
            "limitations": ["Legacy receipt does not establish validated analysis or supported interpretation."]}


def inspect_run(run_root: Path, *, verify: bool = False) -> dict[str, Any]:
    stages = []
    for path in sorted((run_root / "state/stages").glob("*.json")):
        receipt = json.loads(path.read_text(encoding="utf-8"))
        issues, unbound = [], []
        records = []
        for group in ("inputs", "outputs"):
            values = receipt.get(group, [])
            records.extend(values.values() if isinstance(values, dict) else values)
        if isinstance(receipt.get("source_manifest"), dict):
            records.append(receipt["source_manifest"])
        for record in records:
            if not isinstance(record, dict) or "path" not in record:
                unbound.append(str(record))
                continue
            artifact = Path(record["path"])
            if not artifact.is_absolute():
                artifact = run_root / artifact
            if not artifact.is_file():
                issues.append({"path": str(artifact), "problem": "missing"})
            elif not record.get("sha256"):
                unbound.append(str(artifact))
            elif verify and sha256_file(artifact) != record["sha256"]:
                issues.append({"path": str(artifact), "problem": "changed"})
        if issues:
            integrity = "stale_or_missing"
        elif not receipt.get("inputs") or not receipt.get("outputs") or not records or unbound:
            integrity = "not_fully_bound"
        else:
            integrity = "verified" if verify else "not_verified"
        stages.append({"stage": receipt.get("stage", path.stem), "route": receipt.get("route", "legacy_managed"),
                       **_assessment(receipt), "integrity": integrity,
                       "artifact_issues": issues, "unbound_artifacts": unbound, "receipt": str(path)})
    return {"schema_version": "scoop.run-status.v1", "run_root": str(run_root.resolve()),
            "verification_scope": "declared_artifact_integrity_only_not_scientific_truth", "stages": stages,
            "stale_or_missing_stages": [s["stage"] for s in stages if s["artifact_issues"]],
            "unassessed_or_limited_stages": [s["stage"] for s in stages if s["validation"]["status"] not in {"passed", "not_applicable"}]}
