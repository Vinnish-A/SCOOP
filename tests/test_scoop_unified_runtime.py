from __future__ import annotations

import subprocess
import sys
import tomllib
from pathlib import Path

from packaging.requirements import Requirement

ROOT = Path(__file__).resolve().parents[1]


def _dependency_names(values: list[str]) -> set[str]:
    return {Requirement(value).name.lower() for value in values}


def test_pyproject_has_one_managed_runtime_extra() -> None:
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    extras = project["project"]["optional-dependencies"]
    assert set(extras) == {"runtime"}
    assert {
        "anndataoom",
        "cellphonedb",
        "cupy-cuda12x",
        "fastccc",
        "liana",
        "omicverse",
        "scarches",
        "scvi-tools",
        "setuptools",
        "spatialdata",
        "squidpy",
    } <= _dependency_names(extras["runtime"])
    assert "dev" in project["dependency-groups"]


def test_legacy_environment_sources_are_absent() -> None:
    assert not (ROOT / "environment.yml").exists()
    assert not (ROOT / "environment_omicverse.yml").exists()
    assert not (ROOT / "third_party" / "omicverse").exists()
    assert "third_party/omicverse" not in (ROOT / ".gitmodules").read_text(encoding="utf-8")


def test_installer_exposes_only_unified_runtime() -> None:
    completed = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "install_scoop.py"), "--help"],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    assert "--profile" not in completed.stdout
    assert "--home" in completed.stdout
