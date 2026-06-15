from __future__ import annotations

from pathlib import Path
from typing import Any
import json

from .storage import sha256_file, write_json


def build_manifest(run_root: str | Path) -> dict[str, Any]:
    run_root = Path(run_root)
    files = []
    for path in sorted(run_root.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            rel = path.relative_to(run_root)
            files.append({
                "path": str(rel),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            })
    return {"run_root": str(run_root), "n_files": len(files), "files": files}


def write_manifest(run_root: str | Path) -> Path:
    manifest = build_manifest(run_root)
    return write_json(manifest, Path(run_root) / "artifacts" / "manifest.json")
