from __future__ import annotations

import zipfile
from pathlib import Path

import pytest

from scoop.fastccc_cli import _DATABASE_FILES, database_directory


def _write_database(root: Path) -> None:
    root.mkdir(parents=True)
    for name in _DATABASE_FILES:
        (root / name).write_text("placeholder\n", encoding="utf-8")


def test_database_directory_accepts_directory(tmp_path: Path) -> None:
    database = tmp_path / "database"
    _write_database(database)
    with database_directory(database) as resolved:
        assert resolved == database.resolve()


def test_database_directory_finds_database_in_zip(tmp_path: Path) -> None:
    database = tmp_path / "source" / "nested"
    _write_database(database)
    archive = tmp_path / "database.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        for path in database.iterdir():
            handle.write(path, f"nested/{path.name}")
    with database_directory(archive) as resolved:
        assert {path.name for path in resolved.iterdir()} == set(_DATABASE_FILES)


def test_database_directory_rejects_zip_traversal(tmp_path: Path) -> None:
    archive = tmp_path / "unsafe.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("../escape.csv", "no")
    with pytest.raises(ValueError, match="unsafe path"):
        with database_directory(archive):
            pass
