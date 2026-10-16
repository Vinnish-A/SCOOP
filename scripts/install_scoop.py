#!/usr/bin/env python3
"""Build, verify, and atomically install SCOOP under SCOOP_HOME."""

from __future__ import annotations

import argparse
import configparser
import hashlib
import json
import os
import platform
import shutil
import subprocess
import tempfile
import uuid
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PYTHON_VERSION = "3.11"
RUNTIME_NAME = "runtime"
SOURCE_DATE_EPOCH = "315532800"
REQUIRED_CAPABILITIES = ("core", "atlas", "spatial", "ccc", "omicverse", "oom", "parallel_leiden", "programmes")
RUNTIME_IMPORTS = (
    "networkit",
    "rctd",
    "cnmf",
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


def _python_identity() -> tuple[str, dict[str, str]]:
    executable = subprocess.check_output(["uv", "python", "find", PYTHON_VERSION], cwd=ROOT, text=True).strip()
    code = (
        "import json,platform,sysconfig; "
        "print(json.dumps(dict(implementation=platform.python_implementation(), "
        "version=platform.python_version(), platform=platform.system().lower(), "
        "machine=platform.machine(), abi=sysconfig.get_config_var('SOABI'))))"
    )
    identity = json.loads(subprocess.check_output([executable, "-c", code], text=True))
    return executable, identity


def _dependency_identity(requirements: Path, python_identity: dict[str, str]) -> str:
    digest = hashlib.sha256(requirements.read_bytes())
    digest.update(json.dumps(python_identity, sort_keys=True).encode())
    return digest.hexdigest()[:20]


def _completion_record(dependency_id: str, python_identity: dict[str, str]) -> dict[str, object]:
    return {
        "dependency_id": dependency_id,
        "runtime": RUNTIME_NAME,
        "python_identity": python_identity,
        "smoke_imports": list(_runtime_imports()),
    }


def _write_application(staging: Path, wheel: Path) -> None:
    # The SCOOP application belongs to its immutable release, not to the
    # reusable dependency environment. Rollback therefore restores code too.
    with zipfile.ZipFile(wheel) as archive:
        archive.extractall(staging / "lib")
    parser = configparser.ConfigParser()
    entry_file = next((staging / "lib").glob("*.dist-info/entry_points.txt"))
    parser.read(entry_file)
    entry_points = dict(parser["console_scripts"])
    bootstrap = (
        """from pathlib import Path
import importlib
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / 'lib'))
ENTRY_POINTS = %r
if __name__ == '__main__':
    entry = sys.argv.pop(1)
    module, function = ENTRY_POINTS[entry].split(':')
    raise SystemExit(getattr(importlib.import_module(module), function)())
"""
        % entry_points
    )
    (staging / "bootstrap.py").write_text(bootstrap)
    bin_dir = staging / "bin"
    bin_dir.mkdir()
    for command in entry_points:
        launcher = bin_dir / command
        launcher.write_text(
            "#!/bin/sh\nset -eu\n"
            'SCOOP_RELEASE=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd -P)\n'
            'SCOOP_ENV=$(cat "$SCOOP_RELEASE/environment.path")\n'
            'export PATH="$SCOOP_RELEASE/bin:$SCOOP_ENV/bin:$PATH"\n'
            'exec "$SCOOP_ENV/bin/python" -I "$SCOOP_RELEASE/bootstrap.py" ' + command + ' "$@"\n'
        )
        launcher.chmod(0o755)
    # First-class API interpreter: share the exact dependency Python while
    # pinning SCOOP imports to this release. No second venv or pip mutation.
    python_launcher = bin_dir / "python"
    python_launcher.write_text(
        "#!/bin/sh\n"
        "set -eu\n"
        'SCOOP_RELEASE=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd -P)\n'
        'SCOOP_ENV=$(cat "$SCOOP_RELEASE/environment.path")\n'
        'export PYTHONPATH="$SCOOP_RELEASE/lib${PYTHONPATH:+:$PYTHONPATH}"\n'
        'export PATH="$SCOOP_RELEASE/bin:$SCOOP_ENV/bin:$PATH"\n'
        'exec "$SCOOP_ENV/bin/python" "$@"\n'
    )
    python_launcher.chmod(0o755)


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
    digest = hashlib.sha256(Path(__file__).read_bytes())
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


def _release_ready(path: Path) -> bool:
    return path.is_dir() and (path / "manifest.json").is_file() and (path / "environment.path").is_file()


def _select_previous_release(home: Path, current: Path, candidate: Path) -> str | None:
    if not current.is_symlink():
        return None
    active = current.resolve()
    if active != candidate:
        return str(active)
    activation = home / "state" / "activation.json"
    if activation.is_file():
        try:
            recorded = json.loads(activation.read_text(encoding="utf-8")).get("previous_release")
        except (OSError, json.JSONDecodeError):
            recorded = None
        if recorded:
            previous = Path(recorded).resolve()
            if previous != candidate and _release_ready(previous):
                return str(previous)
    alternatives = sorted(
        (path for path in (home / "releases").iterdir() if path != candidate and _release_ready(path)),
        key=lambda path: path.stat().st_mtime_ns,
        reverse=True,
    )
    return str(alternatives[0].resolve()) if alternatives else None


def _write_launcher(home: Path) -> None:
    entry_points = sorted(path.name for path in (home / "current" / "bin").iterdir())
    for command in entry_points:
        launcher = home / "bin" / command
        launcher.parent.mkdir(parents=True, exist_ok=True)
        content = """#!/bin/sh
set -eu
if [ -n "${SCOOP_HOME:-}" ]; then
    SCOOP_ROOT=$SCOOP_HOME
else
    SCOOP_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd -P)
fi
export SCOOP_HOME="$SCOOP_ROOT"
if [ -x "$SCOOP_ROOT/current/bin/%s" ]; then
    exec "$SCOOP_ROOT/current/bin/%s" "$@"
fi
exec "$(cat "$SCOOP_ROOT/current/environment.path")/bin/%s" "$@"
""" % (command, command, command)
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
# Verify the AnnData/SciPy storage contract, not only successful imports.
import tempfile
from pathlib import Path
import anndata as ad
import numpy as np
from scipy import sparse
with tempfile.TemporaryDirectory(prefix="scoop-storage-smoke-") as folder:
    path = Path(folder) / "counts.h5ad"
    ad.AnnData(sparse.csr_matrix([[0, 1], [2, 0]])).write_h5ad(path)
    backed = ad.read_h5ad(path, backed="r")
    try:
        np.testing.assert_array_equal(backed[[1, 0, 1], :].X.toarray(), [[2, 0], [0, 1], [2, 0]])
    finally:
        backed.file.close()
from cnmf.cnmf import cNMF
import networkit as nk
if not hasattr(nk.community, "ParallelLeiden"):
    raise RuntimeError("NetworKit ParallelLeiden API is unavailable")
import rctd
if not all(hasattr(rctd, key) for key in ("Reference", "RCTD", "run_full_mode", "run_multi_mode")):
    raise RuntimeError("RCTD public APIs are unavailable")
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


def install(home: Path, *, register_skill: bool = True) -> Path:
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
    if register_skill and codex_link.exists() and not codex_link.is_symlink():
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
            "--no-header",
            "--no-annotate",
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
        python_executable, python_identity = _python_identity()
        dependency_id = _dependency_identity(requirements, python_identity)
        content_id = _content_id(skill_source, wheel, lock, requirements)
        release = home / "releases" / content_id
        environment = home / "envs" / f"{RUNTIME_NAME}-{dependency_id}"
        current = home / "current"

        if not release.exists():
            staging = home / f".release-{content_id}-{uuid.uuid4().hex}"
            staging.mkdir()
            shutil.copytree(skill_source, staging / "skill")
            shutil.copy2(wheel, staging / wheel.name)
            shutil.copy2(lock, staging / "uv.lock")
            shutil.copy2(requirements, staging / requirements.name)
            _write_application(staging, wheel)
            shutil.copy2(Path(__file__), staging / "install_scoop.py")
            (staging / "environment.path").write_text(str(environment.resolve()) + "\n")
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
                "schema_version": "scoop.release.v3",
                "dependency_id": dependency_id,
                "python_identity": python_identity,
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

        previous = _select_previous_release(home, current, release)

        completion_path = environment / ".complete.json"
        expected_completion = _completion_record(dependency_id, python_identity)
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
            _run(["uv", "venv", "--python", python_executable, str(environment)])
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
            _smoke_runtime(environment / "bin" / "python")
            completion_path.write_text(json.dumps(expected_completion) + "\n", encoding="utf-8")

        process_environment = dict(os.environ)
        process_environment["SCOOP_HOME"] = str(home)
        _run([str(release / "bin" / "scoop"), "--version"], environment=process_environment)
        doctor_process = subprocess.run(
            [str(release / "bin" / "scoop"), "doctor", "--json"],
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
        if register_skill:
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
    parser.add_argument(
        "--no-skill-link",
        action="store_true",
        help="Install an isolated runtime without changing the global Skill link",
    )
    args = parser.parse_args()
    print(install(Path(args.home).expanduser().resolve(), register_skill=not args.no_skill_link))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
