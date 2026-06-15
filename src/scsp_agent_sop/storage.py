from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Mapping
import hashlib
import json
import os

import pandas as pd


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def write_table(df: pd.DataFrame, path: str | Path, index: bool = False) -> Path:
    path = Path(path)
    ensure_dir(path.parent)
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        df.to_parquet(path, index=index)
    elif suffix in {".tsv", ".txt"}:
        df.to_csv(path, sep="\t", index=index)
    elif suffix == ".csv":
        df.to_csv(path, index=index)
    else:
        raise ValueError(f"Unsupported table format: {path}")
    return path


def write_json(obj: Any, path: str | Path) -> Path:
    path = Path(path)
    ensure_dir(path.parent)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return path


def append_jsonl(obj: Mapping[str, Any], path: str | Path) -> Path:
    path = Path(path)
    ensure_dir(path.parent)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(obj, ensure_ascii=False, sort_keys=True) + "\n")
    return path


def init_file_registry(adata, run_id: str | None = None) -> None:
    adata.uns.setdefault("file_registry", {})
    if run_id is not None:
        adata.uns["file_registry"].setdefault("run_id", run_id)
    adata.uns["file_registry"].setdefault("tables", {})
    adata.uns["file_registry"].setdefault("figures", {})
    adata.uns["file_registry"].setdefault("reports", {})
    adata.uns["file_registry"].setdefault("artifacts", {})


def register_file(
    adata,
    *,
    key: str,
    path: str | Path,
    category: str = "tables",
    schema: str | None = None,
    extra: Mapping[str, Any] | None = None,
) -> None:
    path = Path(path)
    init_file_registry(adata)
    rec: dict[str, Any] = {
        "path": os.fspath(path),
        "format": path.suffix.lstrip("."),
    }
    if path.exists() and path.is_file():
        rec["sha256"] = sha256_file(path)
        rec["bytes"] = path.stat().st_size
    if schema:
        rec["schema"] = schema
    if extra:
        rec.update(dict(extra))
    adata.uns["file_registry"].setdefault(category, {})[key] = rec


def flatten_keep_schema(schema: Mapping[str, Any], key: str) -> list[str]:
    item = schema.get(key, [])
    if isinstance(item, list):
        return list(item)
    if isinstance(item, dict):
        out: list[str] = []
        for values in item.values():
            if isinstance(values, list):
                out.extend(values)
        return out
    return []


def prune_mapping(mapping, keep: Iterable[str]) -> None:
    keep_set = set(keep)
    for k in list(mapping.keys()):
        if k not in keep_set:
            del mapping[k]


def prune_h5ad(adata, schema: Mapping[str, Any]) -> None:
    """Prune an AnnData object according to h5ad_schema.yaml.

    This function is deliberately conservative: missing keys are ignored,
    existing allowed keys are retained, and everything else is removed from
    obs/var/obsm/obsp/layers/uns.
    """
    obs_keep = [c for c in flatten_keep_schema(schema, "obs_keep") if c in adata.obs]
    var_keep = [c for c in schema.get("var_keep", []) if c in adata.var]
    adata.obs.drop(columns=[c for c in adata.obs.columns if c not in obs_keep], inplace=True)
    adata.var.drop(columns=[c for c in adata.var.columns if c not in var_keep], inplace=True)
    prune_mapping(adata.layers, schema.get("layers_keep", []))
    prune_mapping(adata.obsm, schema.get("obsm_keep", []))
    prune_mapping(adata.obsp, schema.get("obsp_keep", []))
    uns_keep = set(schema.get("uns_keep", []))
    for key in list(adata.uns.keys()):
        if key not in uns_keep:
            del adata.uns[key]
