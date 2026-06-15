from __future__ import annotations

from typing import Any


def require_omicverse():
    try:
        import omicverse as ov
    except Exception as exc:
        raise ImportError("OmicVerse is not installed. Install `omicverse` or disable this optional facility.") from exc
    return ov


def read_h5ad(path: str):
    ov = require_omicverse()
    return ov.io.read_h5ad(path)


def save(adata, path: str):
    ov = require_omicverse()
    return ov.io.save(adata, path)


def to_gpu(adata):
    ov = require_omicverse()
    return ov.pp.anndata_to_GPU(adata)


def to_cpu(adata):
    ov = require_omicverse()
    return ov.pp.anndata_to_CPU(adata)


def generate_report(adata, output_path: str, species: str = "human", sample_key: str | None = None):
    ov = require_omicverse()
    return ov.single.generate_scRNA_report(adata, output_path=output_path, species=species, sample_key=sample_key)
