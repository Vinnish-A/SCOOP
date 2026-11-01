from __future__ import annotations

from pathlib import Path
from typing import Any
import json

from .storage import sha256_file, write_json


def build_manifest(run_root: str | Path, *, paths=None, metadata=None) -> dict[str, Any]:
    run_root = Path(run_root).resolve()
    files = []
    selected = sorted(run_root.rglob("*")) if paths is None else sorted({Path(p).resolve() for p in paths})
    missing = []
    for path in selected:
        if not path.is_file():
            missing.append(str(path))
            continue
        if path.is_file() and path.name != "manifest.json":
            rel = path.relative_to(run_root.resolve()) if path.is_relative_to(run_root.resolve()) else path
            files.append({
                "path": str(rel),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            })
    return {"run_root": str(run_root), "n_files": len(files), "files": files,
            "scope": "registered_delivery_artifacts" if paths is not None else "run_directory_inventory",
            "missing_registered_artifacts": missing, **(metadata or {})}


def write_manifest(run_root: str | Path, *, paths=None, metadata=None) -> Path:
    manifest = build_manifest(run_root, paths=paths, metadata=metadata)
    return write_json(manifest, Path(run_root) / "artifacts" / "manifest.json")
