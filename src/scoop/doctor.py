"""Capability-aware diagnostics for the managed SCOOP runtime."""

from __future__ import annotations

import json
import os
import platform
import shutil
import sys
from importlib import metadata, util
from pathlib import Path
from typing import Any

from .home import ScoopHomeError, resolve_home

MODULE_NAMES = {
    "pyyaml": "yaml",
    "scvi-tools": "scvi",
}


def _package(name: str, *, required: bool) -> dict[str, Any]:
    available = util.find_spec(MODULE_NAMES.get(name, name.replace("-", "_"))) is not None
    version: str | None = None
    if available:
        try:
            version = metadata.version(name)
        except metadata.PackageNotFoundError:
            version = "unknown"
    return {"available": available, "required": required, "version": version}


def _command(name: str) -> dict[str, Any]:
    path = shutil.which(name)
    return {"available": path is not None, "path": path}


def health() -> dict[str, Any]:
    try:
        home = resolve_home(require_exists=False)
        home_error = None
    except ScoopHomeError as exc:
        home = Path.home() / ".scoop"
        home_error = str(exc)

    packages = {
        name: _package(name, required=True) for name in ("anndata", "scanpy", "numpy", "pandas", "scipy", "pyyaml")
    }
    optional = {
        name: _package(name, required=False)
        for name in ("scvi-tools", "scarches", "omicverse", "squidpy", "spatialdata", "liana", "cellphonedb")
    }
    core_ok = all(value["available"] for value in packages.values())
    atlas_ok = optional["scvi-tools"]["available"] and optional["scarches"]["available"]
    try:
        import torch

        torch_check = {
            "available": True,
            "version": torch.__version__,
            "cuda_available": bool(torch.cuda.is_available()),
            "cuda_version": torch.version.cuda,
        }
    except Exception as exc:  # pragma: no cover - environment dependent
        torch_check = {"available": False, "error": str(exc)}

    return {
        "schema_version": "scoop.doctor.v1",
        "ok": core_ok and home_error is None,
        "home": str(home),
        "home_error": home_error,
        "runtime": {
            "python": sys.version.split()[0],
            "executable": sys.executable,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cwd": os.getcwd(),
        },
        "capabilities": {
            "core": {"ok": core_ok, "packages": packages},
            "atlas": {"ok": atlas_ok, "packages": {k: optional[k] for k in ("scvi-tools", "scarches")}},
            "spatial": {"ok": optional["squidpy"]["available"] and optional["spatialdata"]["available"]},
            "ccc": {"ok": optional["liana"]["available"] or optional["cellphonedb"]["available"]},
            "omicverse": {"ok": optional["omicverse"]["available"]},
            "r-reference": {"ok": _command("Rscript")["available"], "Rscript": _command("Rscript")},
            "gpu": torch_check,
        },
        "commands": {name: _command(name) for name in ("uv", "Rscript", "git")},
    }


def emit(*, json_output: bool) -> int:
    result = health()
    if json_output:
        print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    else:
        print(f"SCOOP Home: {result['home']}")
        for name, capability in result["capabilities"].items():
            print(f"{name}: {'ok' if capability.get('ok') else 'unavailable'}")
    return 0 if result["ok"] else 1
