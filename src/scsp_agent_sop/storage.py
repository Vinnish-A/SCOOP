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


def _write_atomic(path: str | Path, write, suffix: str) -> Path:
    """A serialization failure must not overwrite a previously valid file."""
    from scoop.io import default_file_mode

    path = Path(path)
    ensure_dir(path.parent)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.stem}-", suffix=suffix, dir=path.parent)
    os.close(fd)
    try:
        write(temporary)
        os.chmod(temporary, default_file_mode())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return path


def write_h5ad_atomic(adata, path: str | Path) -> Path:
    return _write_atomic(path, adata.write_h5ad, ".h5ad")


def _obs_names_sha256(obs_names) -> str:
    return hashlib.sha256("\n".join(map(str, obs_names)).encode("utf-8")).hexdigest()


def write_obs_sidecar(adata, columns: Iterable[str], path: str | Path, *, base: str | Path,
                      key: str | None = None) -> Path:
    """Write obs_name plus ``columns`` to parquet instead of a full H5AD copy (DC1, C1-08, S032).

    ``columns`` are the obs columns the stage declares as its output, whether or
    not a value differs from ``base``, the H5AD or sidecar the stage read. The
    table is registered under ``key`` (default: the file stem). A JSON receipt
    beside it binds the table and its base by sha256 and records the obs_name
    order and the file registry an H5AD copy would have carried in uns;
    read_adata overlays the table onto that base. Both files are replaced
    atomically. Returns the receipt path.
    """
    from scoop.io import atomic_json

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
    table = table.reset_index(names="obs_name")
    _write_atomic(path, lambda temporary: table.to_parquet(temporary, index=False), ".parquet")
    base = Path(base).resolve()
    register_file(adata, key=key, path=path, schema="obs_sidecar.v2",
                  extra={"columns": columns, "n_obs": int(adata.n_obs),
                         "base": {"path": str(base), "sha256": sha256_file(base),
                                  "obs_names_sha256": _obs_names_sha256(adata.obs_names)}})
    registry = adata.uns["file_registry"]
    receipt = {"schema_version": "scoop.obs_sidecar_receipt.v2", "sidecar": registry["tables"][key], "file_registry": registry}
    receipt = json.loads(json.dumps(receipt, default=lambda v: v.tolist() if hasattr(v, "tolist") else str(v)))
    return atomic_json(path.with_suffix(".json"), receipt)


def read_adata(path: str | Path):
    """Read an H5AD, or an obs sidecar overlaid onto the base it was written against.

    The overlay is the stage's obs and file registry on the base's matrices and
    other slots. It refuses a sidecar whose table, base file or base obs_name
    order differs from what its receipt recorded.
    """
    import anndata as ad
    import numpy as np

    path = Path(path)
    if path.suffix != ".parquet":
        return ad.read_h5ad(path)
    receipt_path = path.with_suffix(".json")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    sidecar = receipt["sidecar"]
    if "base" not in sidecar:
        raise ValueError(f"{receipt_path} records no base input; rerun the stage that wrote {path.name}.")
    if sha256_file(path) != sidecar["sha256"]:
        raise ValueError(f"{path} differs from the table its receipt {receipt_path} records.")
    base = sidecar["base"]
    if sha256_file(base["path"]) != base["sha256"]:
        raise ValueError(f"{base['path']} changed after the sidecar {path} was written against it; "
                         "rerun the stage that wrote the sidecar.")
    adata = read_adata(base["path"])
    if _obs_names_sha256(adata.obs_names) != base["obs_names_sha256"]:
        raise ValueError(f"{path} rows do not follow the observation order of its base {base['path']}.")
    table = pd.read_parquet(path, columns=sidecar["columns"]).set_axis(adata.obs_names)
    for column in sidecar["columns"]:
        values = table[column]
        if isinstance(values.dtype, pd.StringDtype):
            # anndata does not write nullable string arrays; restore the object column
            values = values.astype(object).where(values.notna(), np.nan)
        adata.obs[column] = values
    adata.uns["file_registry"] = receipt["file_registry"]
    return adata


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
    obs/var/obsm/obsp/layers/uns. Only the schema names retained keys; uns
    graph, UMAP/PCA and log1p descriptors follow their retained matrices.
    """
    obs_keep = set(flatten_keep_schema(schema, "obs_keep"))
    var_keep = set(schema.get("var_keep", []))
    adata.obs.drop(columns=[c for c in adata.obs.columns if c not in obs_keep], inplace=True)
    adata.var.drop(columns=[c for c in adata.var.columns if c not in var_keep], inplace=True)
    prune_mapping(adata.layers, schema.get("layers_keep", []))
    prune_mapping(adata.obsm, schema.get("obsm_keep", []))
    prune_mapping(adata.obsp, schema.get("obsp_keep", []))
    uns_keep = set(schema.get("uns_keep", []))
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
