#!/usr/bin/env python3
"""Build, verify, and atomically install SCOOP under SCOOP_HOME."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import tempfile
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PYTHON_VERSION = "3.11"
RUNTIME_NAME = "runtime"
SOURCE_DATE_EPOCH = "315532800"
REQUIRED_CAPABILITIES = ("core", "atlas", "spatial", "ccc", "omicverse", "oom")
RUNTIME_IMPORTS = (
    "anndataoom",
    "cellphonedb",
    "fastccc",
    "liana",
    "omicverse",
    "scarches",
    "scvi",
    "spatialdata",
    "squidpy",
)


def _runtime_imports() -> tuple[str, ...]:
    imports: tuple[str, ...] = RUNTIME_IMPORTS
    if platform.system() == "Linux" and platform.machine() == "x86_64":
        imports += ("cupy",)
    return imports


def _completion_record(content_id: str) -> dict[str, object]:
    return {
        "content_id": content_id,
        "runtime": RUNTIME_NAME,
        "python": PYTHON_VERSION,
        "platform": platform.system().lower(),
        "machine": platform.machine(),
        "smoke_imports": list(_runtime_imports()),
    }


def _run(arguments: list[str], *, cwd: Path = ROOT, environment: dict[str, str] | None = None) -> None:
    subprocess.run(arguments, cwd=cwd, check=True, env=environment)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _files(root: Path) -> list[Path]:
    return sorted(path for path in root.rglob("*") if path.is_file() and "__pycache__" not in path.parts)


def _content_id(skill: Path, wheel: Path, lock: Path, requirements: Path) -> str:
    digest = hashlib.sha256()
    for base, paths in (
        (skill, _files(skill)),
        (wheel.parent, [wheel]),
        (lock.parent, [lock]),
        (requirements.parent, [requirements]),
    ):
        for path in paths:
            digest.update(path.relative_to(base).as_posix().encode())
            digest.update(_sha256(path).encode())
    return digest.hexdigest()[:20]


def _atomic_symlink(target: Path, link: Path) -> None:
    link.parent.mkdir(parents=True, exist_ok=True)
    temporary = link.parent / f".{link.name}.{uuid.uuid4().hex}.tmp"
    temporary.symlink_to(target)
    os.replace(temporary, link)


def _write_launcher(home: Path) -> None:
    launcher = home / "bin" / "scoop"
    launcher.parent.mkdir(parents=True, exist_ok=True)
    content = """#!/bin/sh
set -eu
SCOOP_ROOT=${SCOOP_HOME:-\"$HOME/.scoop\"}
exec \"$(cat \"$SCOOP_ROOT/current/environment.path\")/bin/scoop\" \"$@\"
"""
    temporary = launcher.with_name(f".{launcher.name}.{uuid.uuid4().hex}.tmp")
    temporary.write_text(content, encoding="utf-8")
    temporary.chmod(0o755)
    os.replace(temporary, launcher)


def _smoke_runtime(python: Path) -> None:
    code = """
import importlib
modules = %r
for name in modules:
    importlib.import_module(name)
import omicverse as ov
if not hasattr(ov, "pp") or not hasattr(ov, "single"):
    raise RuntimeError("OmicVerse public APIs pp/single are unavailable")
""" % (_runtime_imports(),)
    completed = subprocess.run(
        [str(python), "-c", code],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if completed.returncode:
        detail = (completed.stderr or completed.stdout)[-4000:]
        raise RuntimeError(f"candidate unified runtime import smoke failed:\n{detail}")


def install(home: Path) -> Path:
    if not home.is_absolute():
        raise ValueError("installation home must be absolute")
    home.mkdir(parents=True, exist_ok=True)
    (home / "releases").mkdir(exist_ok=True)
    (home / "envs").mkdir(exist_ok=True)
    skill_source = ROOT / "skill" / "scoop"
    if not (skill_source / "SKILL.md").is_file():
        raise FileNotFoundError("skill/scoop/SKILL.md")
    lock = ROOT / "uv.lock"
    if not lock.is_file():
        raise FileNotFoundError("uv.lock; run uv lock before installation")

    codex_link = Path.home() / ".codex" / "skills" / "scoop"
    if codex_link.exists() and not codex_link.is_symlink():
        raise RuntimeError(f"refusing to replace non-symlink Codex Skill path: {codex_link}")

    with tempfile.TemporaryDirectory(prefix="scoop-build-") as temporary_name:
        temporary = Path(temporary_name)
        dist = temporary / "dist"
        build_environment = dict(os.environ)
        build_environment["SOURCE_DATE_EPOCH"] = SOURCE_DATE_EPOCH
        _run(["uv", "build", "--quiet", "--wheel", "--out-dir", str(dist)], environment=build_environment)
        wheel = next(dist.glob("*.whl"))
        requirements = temporary / "requirements-runtime.lock"
        export = [
            "uv",
            "export",
            "--frozen",
            "--no-dev",
            "--no-emit-project",
            "--extra",
            "runtime",
        ]
        export_process = subprocess.run(
            export,
            cwd=ROOT,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
        )
        requirements.write_text(export_process.stdout, encoding="utf-8")
        content_id = _content_id(skill_source, wheel, lock, requirements)
        release = home / "releases" / content_id
        environment = home / "envs" / f"{RUNTIME_NAME}-{content_id}"
        current = home / "current"
        previous = str(current.resolve()) if current.is_symlink() else None

        if not release.exists():
            staging = home / f".release-{content_id}-{uuid.uuid4().hex}"
            staging.mkdir()
            shutil.copytree(skill_source, staging / "skill")
            shutil.copy2(wheel, staging / wheel.name)
            shutil.copy2(lock, staging / "uv.lock")
            shutil.copy2(requirements, staging / requirements.name)
            records = []
            for path in _files(staging):
                records.append(
                    {
                        "path": path.relative_to(staging).as_posix(),
                        "sha256": _sha256(path),
                        "bytes": path.stat().st_size,
                    }
                )
            manifest = {
                "schema_version": "scoop.release.v2",
                "content_id": content_id,
                "runtime": RUNTIME_NAME,
                "python": PYTHON_VERSION,
                "platform": platform.system().lower(),
                "machine": platform.machine(),
                "files": records,
            }
            (staging / "manifest.json").write_text(
                json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
            )
            os.replace(staging, release)

        completion_path = environment / ".complete.json"
        expected_completion = _completion_record(content_id)
        if completion_path.is_file():
            try:
                actual_completion = json.loads(completion_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                raise RuntimeError(f"invalid environment completion manifest: {completion_path}") from exc
            if actual_completion != expected_completion:
                raise RuntimeError(f"environment identity mismatch: {environment}")
        else:
            if environment.exists():
                raise RuntimeError(f"incomplete environment requires manual inspection: {environment}")
            _run(["uv", "venv", "--python", PYTHON_VERSION, str(environment)])
            release_requirements = release / "requirements-runtime.lock"
            _run(
                [
                    "uv",
                    "pip",
                    "install",
                    "--python",
                    str(environment / "bin" / "python"),
                    "-r",
                    str(release_requirements),
                ]
            )
            release_wheel = next(release.glob("*.whl"))
            _run(
                [
                    "uv",
                    "pip",
                    "install",
                    "--python",
                    str(environment / "bin" / "python"),
                    "--no-deps",
                    str(release_wheel),
                ]
            )
            _run([str(environment / "bin" / "scoop"), "--version"])
            _smoke_runtime(environment / "bin" / "python")
            completion_path.write_text(json.dumps(expected_completion) + "\n", encoding="utf-8")

        (release / "environment.path").write_text(str(environment.resolve()) + "\n", encoding="utf-8")
        process_environment = dict(os.environ)
        process_environment["SCOOP_HOME"] = str(home)
        doctor_process = subprocess.run(
            [str(environment / "bin" / "scoop"), "doctor", "--json"],
            cwd=home,
            check=True,
            env=process_environment,
            text=True,
            stdout=subprocess.PIPE,
        )
        doctor = json.loads(doctor_process.stdout)
        unavailable = [name for name in REQUIRED_CAPABILITIES if not doctor["capabilities"][name]["ok"]]
        if unavailable:
            raise RuntimeError(f"candidate unified runtime lacks capabilities: {', '.join(unavailable)}")
        _atomic_symlink(release, current)
        _atomic_symlink(release / "skill", home / "skills" / "scoop" / "current")
        _write_launcher(home)
        codex_link.parent.mkdir(parents=True, exist_ok=True)
        _atomic_symlink(home / "skills" / "scoop" / "current", codex_link)
        activation = home / "state" / "activation.json"
        activation.parent.mkdir(parents=True, exist_ok=True)
        activation_temp = activation.with_name(f".{activation.name}.{uuid.uuid4().hex}.tmp")
        activation_temp.write_text(
            json.dumps({"current_release": str(release), "previous_release": previous}, indent=2) + "\n",
            encoding="utf-8",
        )
        os.replace(activation_temp, activation)
        _run([str(home / "bin" / "scoop"), "doctor", "--json"], cwd=home, environment=process_environment)
        return release


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--home", default=os.environ.get("SCOOP_HOME", str(Path.home() / ".scoop")))
    args = parser.parse_args()
    print(install(Path(args.home).expanduser().resolve()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
