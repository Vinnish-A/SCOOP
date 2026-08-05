"""Stable SCOOP command-line adapter for the FastCCC Python package."""

from __future__ import annotations

import argparse
import shutil
import tempfile
import zipfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

import anndata as ad

_DATABASE_FILES = (
    "complex_composition_table.csv",
    "complex_table.csv",
    "gene_table.csv",
    "interaction_table.csv",
    "protein_table.csv",
)


def _validate_database(path: Path) -> Path:
    missing = [name for name in _DATABASE_FILES if not (path / name).is_file()]
    if missing:
        raise ValueError(f"FastCCC database is missing files under {path}: {', '.join(missing)}")
    return path


def _extract_zip_safely(archive: Path, destination: Path) -> None:
    root = destination.resolve()
    with zipfile.ZipFile(archive) as handle:
        for member in handle.infolist():
            target = (root / member.filename).resolve()
            if not target.is_relative_to(root):
                raise ValueError(f"unsafe path in FastCCC database archive: {member.filename}")
            if member.is_dir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            with handle.open(member) as source, target.open("wb") as sink:
                shutil.copyfileobj(source, sink)


@contextmanager
def database_directory(path: str | Path) -> Iterator[Path]:
    source = Path(path).expanduser().resolve()
    if source.is_dir():
        yield _validate_database(source)
        return
    if not source.is_file():
        raise FileNotFoundError(f"FastCCC database not found: {source}")
    if source.suffix.lower() != ".zip":
        raise ValueError(f"FastCCC database must be a directory or ZIP archive: {source}")
    with tempfile.TemporaryDirectory(prefix="scoop-fastccc-") as temporary:
        root = Path(temporary)
        _extract_zip_safely(source, root)
        candidates = sorted({item.parent for item in root.rglob("interaction_table.csv")})
        valid = [candidate for candidate in candidates if all((candidate / name).is_file() for name in _DATABASE_FILES)]
        if len(valid) != 1:
            raise ValueError(f"expected one FastCCC database in {source}, found {len(valid)}")
        yield valid[0]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run locked FastCCC through SCOOP's stable adapter.")
    parser.add_argument("--input", required=True, help="Annotated H5AD input.")
    parser.add_argument("--groupby", required=True, help="Cell-type column in adata.obs.")
    parser.add_argument("--database", required=True, help="FastCCC/CellPhoneDB database directory or ZIP archive.")
    parser.add_argument("--output", required=True, help="Output directory for FastCCC tables.")
    parser.add_argument("--min-percentile", type=float, default=0.1)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    adata = ad.read_h5ad(args.input)
    if args.groupby not in adata.obs:
        raise KeyError(f"FastCCC groupby column not found in adata.obs: {args.groupby}")
    output = Path(args.output).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    try:
        from fastccc.core import statistical_analysis_method
    except ImportError as exc:
        raise RuntimeError("FastCCC is unavailable; install SCOOP's unified runtime") from exc
    with database_directory(args.database) as database:
        statistical_analysis_method(
            database_file_path=str(database),
            celltype_file_path=None,
            counts_file_path=adata,
            meta_key=args.groupby,
            min_percentile=args.min_percentile,
            save_path=str(output),
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
