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


def test_dependency_identity_ignores_release_code_and_tracks_python(tmp_path):
    requirements = tmp_path / "requirements-runtime.lock"
    requirements.write_text("numpy==2.1.0 --hash=sha256:example\n")
    identity = {"implementation": "CPython", "version": "3.11.12", "machine": "x86_64", "platform": "linux"}
    original = INSTALLER._dependency_identity(requirements, identity)
    # Skill/application files do not participate in dependency identity.
    (tmp_path / "SKILL.md").write_text("changed instructions")
    assert INSTALLER._dependency_identity(requirements, identity) == original
    assert INSTALLER._dependency_identity(requirements, dict(identity, version="3.11.13")) != original
    requirements.write_text("numpy==2.2.0 --hash=sha256:different\n")
    assert INSTALLER._dependency_identity(requirements, identity) != original


def test_release_code_and_rollback_share_dependencies_without_mutating_them(tmp_path):
    import subprocess
    import sys
    import zipfile

    environment = tmp_path / "env"
    (environment / "bin").mkdir(parents=True)
    (environment / "bin/python").symlink_to(sys.executable)
    home = tmp_path / "scoop"
    home.mkdir()
    releases = []
    for name in ("first", "second"):
        release = home / name
        release.mkdir()
        wheel = tmp_path / f"{name}.whl"
        with zipfile.ZipFile(wheel, "w") as archive:
            archive.writestr(
                "scoop_operating_protocol.dist-info/entry_points.txt",
                "[console_scripts]\nscoop = test_application:main\n",
            )
            archive.writestr("test_application.py", f"def main():\n    print({name!r})\n    return 0\n")
        INSTALLER._write_application(release, wheel)
        (release / "environment.path").write_text(str(environment) + "\n")
        releases.append(release)
    INSTALLER._atomic_symlink(releases[0], home / "current")
    INSTALLER._write_launcher(home)
    process_environment = dict(__import__("os").environ, SCOOP_HOME=str(home))

    def launch():
        return subprocess.check_output([str(home / "bin/scoop")], env=process_environment, text=True).strip()

    assert launch() == "first"
    application_path = subprocess.check_output(
        [str(home / "bin/python"), "-c", "import test_application;print(test_application.__file__)"],
        env=process_environment,
        text=True,
    ).strip()
    assert Path(application_path).parent == releases[0] / "lib"
    INSTALLER._atomic_symlink(releases[1], home / "current")
    assert launch() == "second"
    INSTALLER._atomic_symlink(releases[0], home / "current")
    assert launch() == "first"
    assert list((environment / "bin").iterdir()) == [environment / "bin/python"]
