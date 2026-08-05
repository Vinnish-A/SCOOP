"""Hash-bound approval artifacts for the scVI annotation to scPoli gate."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

from .io import atomic_text, sha256_file

SCHEMA = "scoop.atlas-approval.v1"


def _artifact(path: str | Path) -> dict[str, Any]:
    resolved = Path(path).resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    return {"path": str(resolved), "sha256": sha256_file(resolved), "bytes": resolved.stat().st_size}


def approve(
    *,
    gate_path: str | Path,
    run_id: str,
    annotation_version: str,
    approved_by: str,
    input_h5ad: str | Path,
    annotation_report: str | Path,
    label_table: str | Path,
    labeled_indices: str | Path,
    notes: str | None = None,
) -> Path:
    payload = {
        "schema_version": SCHEMA,
        "run_id": run_id,
        "annotation_version": annotation_version,
        "status": "approved",
        "approved_by": approved_by,
        "approved_at": datetime.now(timezone.utc).isoformat(),
        "artifacts": {
            "input_h5ad": _artifact(input_h5ad),
            "annotation_report": _artifact(annotation_report),
            "label_table": _artifact(label_table),
            "labeled_indices": _artifact(labeled_indices),
        },
        "notes": notes,
    }
    return atomic_text(Path(gate_path), yaml.safe_dump(payload, sort_keys=False, allow_unicode=True))


def verify(gate_path: str | Path) -> dict[str, Any]:
    path = Path(gate_path).resolve()
    if not path.is_file():
        return {"ok": False, "errors": [f"gate is missing: {path}"]}
    try:
        loaded = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError as exc:
        return {"ok": False, "errors": [f"invalid gate YAML: {exc}"]}
    if not isinstance(loaded, dict):
        return {"ok": False, "errors": ["approval YAML root must be a mapping"]}
    payload: dict[str, Any] = loaded
    errors: list[str] = []
    if payload.get("schema_version") != SCHEMA:
        errors.append("unsupported approval schema")
    if payload.get("status") != "approved":
        errors.append("approval status is not approved")
    for key in ("approved_by", "approved_at", "run_id", "annotation_version"):
        if not payload.get(key):
            errors.append(f"missing approval field: {key}")
    artifacts_value = payload.get("artifacts")
    artifacts: dict[str, Any] = artifacts_value if isinstance(artifacts_value, dict) else {}
    for name in ("input_h5ad", "annotation_report", "label_table", "labeled_indices"):
        record_value = artifacts.get(name)
        record: dict[str, Any] = record_value if isinstance(record_value, dict) else {}
        artifact_path = Path(str(record.get("path", "")))
        if not artifact_path.is_file():
            errors.append(f"approval artifact missing: {name}")
            continue
        expected = record.get("sha256")
        actual = sha256_file(artifact_path)
        if not expected or expected != actual:
            errors.append(f"approval artifact hash mismatch: {name}")
    return {"ok": not errors, "errors": errors, "gate": str(path), "payload": payload}
