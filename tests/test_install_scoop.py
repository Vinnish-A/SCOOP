from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("install_scoop", ROOT / "scripts" / "install_scoop.py")
assert SPEC is not None and SPEC.loader is not None
INSTALLER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(INSTALLER)


def _ready_release(path: Path) -> None:
    path.mkdir(parents=True)
    (path / "manifest.json").write_text("{}\n", encoding="utf-8")
    (path / "environment.path").write_text("/runtime\n", encoding="utf-8")


def test_idempotent_install_preserves_previous_release(tmp_path: Path) -> None:
    releases = tmp_path / "releases"
    previous = releases / "previous"
    current_release = releases / "current"
    _ready_release(previous)
    _ready_release(current_release)
    (tmp_path / "current").symlink_to(current_release)
    activation = tmp_path / "state" / "activation.json"
    activation.parent.mkdir()
    activation.write_text(json.dumps({"previous_release": str(previous)}) + "\n", encoding="utf-8")

    selected = INSTALLER._select_previous_release(tmp_path, tmp_path / "current", current_release)

    assert selected == str(previous.resolve())


def test_idempotent_install_recovers_latest_alternative(tmp_path: Path) -> None:
    releases = tmp_path / "releases"
    previous = releases / "previous"
    current_release = releases / "current"
    _ready_release(previous)
    _ready_release(current_release)
    (tmp_path / "current").symlink_to(current_release)

    selected = INSTALLER._select_previous_release(tmp_path, tmp_path / "current", current_release)

    assert selected == str(previous.resolve())
