"""Installed release inspection and rollback."""

from __future__ import annotations

import json
import os
import uuid
from pathlib import Path

from .home import resolve_home


def _atomic_symlink(target: Path, link: Path) -> None:
    temporary = link.parent / f".{link.name}.{uuid.uuid4().hex}.tmp"
    temporary.symlink_to(target)
    os.replace(temporary, link)


def status() -> dict[str, object]:
    home = resolve_home()
    current = home / "current"
    if not current.is_symlink():
        return {"ok": False, "home": str(home), "error": "current release link is missing"}
    release = current.resolve()
    manifest_path = release / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
    activation_path = home / "state" / "activation.json"
    activation = json.loads(activation_path.read_text(encoding="utf-8")) if activation_path.is_file() else {}
    return {
        "ok": manifest_path.is_file(),
        "home": str(home),
        "release": str(release),
        "manifest": manifest,
        "activation": activation,
    }


def rollback() -> Path:
    result = status()
    if not result["ok"]:
        raise RuntimeError(str(result.get("error", "active release has no manifest")))
    activation = result["activation"]
    assert isinstance(activation, dict)
    previous = activation.get("previous_release")
    if not previous:
        raise RuntimeError("active release has no rollback target")
    target = Path(str(previous)).resolve()
    if not (target / "manifest.json").is_file():
        raise FileNotFoundError(f"rollback release is unavailable: {target}")
    home = resolve_home()
    current_release = Path(str(result["release"])).resolve()
    _atomic_symlink(target, home / "current")
    _atomic_symlink(target / "skill", home / "skills" / "scoop" / "current")
    activation_path = home / "state" / "activation.json"
    activation_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = activation_path.with_name(f".{activation_path.name}.tmp")
    temporary.write_text(
        json.dumps({"current_release": str(target), "previous_release": str(current_release)}, indent=2) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, activation_path)
    return target
