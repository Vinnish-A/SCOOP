from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Mapping
import hashlib
import json
import os
import tempfile

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


def write_h5ad_atomic(adata, path: str | Path) -> Path:
    """A serialization failure must not overwrite a previously valid object."""
    from scoop.io import default_file_mode

    path = Path(path)
    ensure_dir(path.parent)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.stem}-", suffix=".h5ad", dir=path.parent)
    os.close(fd)
    try:
        adata.write_h5ad(temporary)
        os.chmod(temporary, default_file_mode())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return path


def changed_obs_columns(before: pd.DataFrame, adata) -> list[str]:
    """obs columns that a stage added or rewrote, relative to an earlier copy of obs."""
    return [c for c in adata.obs.columns if c not in before or not adata.obs[c].equals(before[c])]


def write_obs_sidecar(adata, columns: Iterable[str], path: str | Path, *, key: str | None = None) -> Path:
    """Write obs_name plus ``columns`` to parquet instead of a full H5AD copy (DC1, C1-08, S032).

    The table is registered under ``key`` (default: the file stem). A JSON
    receipt beside it records the columns, the obs_name order hash and the
    file registry that an H5AD copy would have carried in uns. Returns the
    receipt path.
    """
    path = Path(path)
    if path.suffix != ".parquet":
        raise ValueError(f"An obs sidecar is a .parquet file: {path}")
    key = key or path.stem
    columns = list(columns)
    table = adata.obs[columns].copy()
    for column in columns:
        # parquet needs one type per column; committed label columns hold mixed objects
        if table[column].dtype == object:
            table[column] = table[column].astype("string")
    write_table(table.reset_index(names="obs_name"), path)
    obs_names_sha256 = hashlib.sha256("\n".join(map(str, adata.obs_names)).encode("utf-8")).hexdigest()
    register_file(adata, key=key, path=path, schema="obs_sidecar.v1",
                  extra={"columns": columns, "n_obs": int(adata.n_obs), "obs_names_sha256": obs_names_sha256})
    registry = adata.uns["file_registry"]
    receipt = {"schema_version": "scoop.obs_sidecar_receipt.v1", "sidecar": registry["tables"][key], "file_registry": registry}
    receipt = json.loads(json.dumps(receipt, default=lambda v: v.tolist() if hasattr(v, "tolist") else str(v)))
    return write_json(receipt, path.with_suffix(".json"))


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
    fingerprint: bool = True,
) -> None:
    path = Path(path)
    init_file_registry(adata)
    rec: dict[str, Any] = {
        "path": os.fspath(path),
        "format": path.suffix.lstrip("."),
    }
    if fingerprint and path.exists() and path.is_file():
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
    obs_keep.extend(c for c in ("annotation_reviewed", "annotation_identity_scope", "broad_annotation_status", "fine_annotation_status") if c in adata.obs)
    var_keep = [c for c in schema.get("var_keep", []) if c in adata.var]
    adata.obs.drop(columns=[c for c in adata.obs.columns if c not in obs_keep], inplace=True)
    adata.var.drop(columns=[c for c in adata.var.columns if c not in var_keep], inplace=True)
    prune_mapping(adata.layers, schema.get("layers_keep", []))
    prune_mapping(adata.obsm, schema.get("obsm_keep", []))
    prune_mapping(adata.obsp, schema.get("obsp_keep", []))
    uns_keep = set(schema.get("uns_keep", []))
    uns_keep.update({"annotation_review", "delivery_scope"})
    # Graph descriptors are executable metadata. Keep descriptors whose matrices
    # survive; old Scanpy-prefixed graphs may have canonical SCOOP copies.
    for name, descriptor in adata.uns.items():
        if not isinstance(descriptor, Mapping) or "connectivities_key" not in descriptor:
            continue
        graph_fields = [field for field in ("connectivities_key", "distances_key") if field in descriptor]
        for field in graph_fields:
            if descriptor[field] not in adata.obsp:
                canonical = f"{field.removesuffix('_key')}_{name.removeprefix('neighbors_')}"
                if canonical in adata.obsp:
                    descriptor[field] = canonical
        if all(descriptor[field] in adata.obsp for field in graph_fields):
            uns_keep.add(name)
    for descriptor, prefix in (("umap", "X_umap"), ("pca", "X_pca")):
        if any(key.startswith(prefix) for key in adata.obsm):
            uns_keep.add(descriptor)
    if "log1p_norm" in adata.layers:
        uns_keep.add("log1p")
    for key in list(adata.uns.keys()):
        if key not in uns_keep:
            del adata.uns[key]
