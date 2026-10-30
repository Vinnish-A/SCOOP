#!/usr/bin/env python3
"""Invoke the managed SCOOP CLI without relying on the caller's PATH."""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


def scoop_home() -> Path:
    raw = os.environ.get("SCOOP_HOME", "").strip()
    home = Path(raw).expanduser() if raw else Path.home() / ".scoop"
    if not home.is_absolute():
        raise SystemExit(f"SCOOP_HOME must be absolute: {home}")
    return home


def main() -> int:
    launcher = scoop_home() / "bin" / "scoop"
    if not launcher.is_file():
        raise SystemExit(f"managed SCOOP launcher is missing: {launcher}")
    return subprocess.run([str(launcher), *sys.argv[1:]], check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
