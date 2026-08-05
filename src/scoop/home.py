"""Resolve the SCOOP installation root without depending on another manager."""

from __future__ import annotations

import os
from pathlib import Path


class ScoopHomeError(ValueError):
    """The configured SCOOP Home is unsafe or unavailable."""


def resolve_home(*, require_exists: bool = True) -> Path:
    raw = os.environ.get("SCOOP_HOME", "").strip()
    candidate = Path(raw).expanduser() if raw else Path.home() / ".scoop"
    if not candidate.is_absolute():
        raise ScoopHomeError(f"SCOOP_HOME must be absolute: {candidate}")
    resolved = candidate.resolve(strict=False)
    if require_exists and not resolved.is_dir():
        raise ScoopHomeError(f"SCOOP Home does not exist: {resolved}")
    return resolved
