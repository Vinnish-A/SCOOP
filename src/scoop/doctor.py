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
    "cupy-cuda12x": "cupy",
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
    extensions = {
        name: _package(name, required=True)
        for name in (
            "anndataoom",
            "cellphonedb",
            "cupy-cuda12x",
            "fastccc",
            "liana",
            "omicverse",
            "scarches",
            "scvi-tools",
            "spatialdata",
            "squidpy",
        )
    }
    core_ok = all(value["available"] for value in packages.values())
    atlas_ok = all(extensions[name]["available"] for name in ("scvi-tools", "scarches"))
    spatial_ok = all(extensions[name]["available"] for name in ("squidpy", "spatialdata"))
    ccc_ok = all(extensions[name]["available"] for name in ("fastccc", "liana", "cellphonedb"))
    omicverse_ok = extensions["omicverse"]["available"]
    oom_ok = omicverse_ok and extensions["anndataoom"]["available"]
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

    torch_check["cupy"] = extensions["cupy-cuda12x"]
    torch_check["ok"] = bool(
        torch_check["available"] and torch_check.get("cuda_available") and extensions["cupy-cuda12x"]["available"]
    )
    unified_ok = core_ok and atlas_ok and spatial_ok and ccc_ok and omicverse_ok and oom_ok

    return {
        "schema_version": "scoop.doctor.v1",
        "ok": unified_ok and home_error is None,
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
            "atlas": {"ok": atlas_ok, "packages": {k: extensions[k] for k in ("scvi-tools", "scarches")}},
            "spatial": {"ok": spatial_ok, "packages": {k: extensions[k] for k in ("squidpy", "spatialdata")}},
            "ccc": {
                "ok": ccc_ok,
                "packages": {k: extensions[k] for k in ("fastccc", "liana", "cellphonedb")},
            },
            "omicverse": {"ok": omicverse_ok, "packages": {"omicverse": extensions["omicverse"]}},
            "oom": {
                "ok": oom_ok,
                "packages": {"anndataoom": extensions["anndataoom"]},
            },
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
