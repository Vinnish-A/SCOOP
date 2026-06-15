from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
import scipy.sparse as sp


def read_counts(path: Path, *, layer: str | None = None, transpose: bool = False) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".h5ad":
        obj = ad.read_h5ad(path)
        matrix = obj.layers[layer] if layer else (obj.layers["counts"] if "counts" in obj.layers else obj.X)
        if sp.issparse(matrix):
            matrix = matrix.toarray()
        frame = pd.DataFrame(np.asarray(matrix), index=obj.obs_names.astype(str), columns=obj.var_names.astype(str))
        if not transpose:
            frame = frame.T
        return frame
    if suffix in {".tsv", ".txt"}:
        frame = pd.read_csv(path, sep="\t", index_col=0)
    elif suffix == ".csv":
        frame = pd.read_csv(path, index_col=0)
    elif suffix in {".parquet", ".pq"}:
        frame = pd.read_parquet(path)
        if frame.columns[0] not in frame.select_dtypes(include="number").columns:
            frame = frame.set_index(frame.columns[0])
    else:
        raise ValueError(f"unsupported counts format: {path}")
    if transpose:
        frame = frame.T
    if not all(is_numeric_dtype(dtype) for dtype in frame.dtypes):
        frame = frame.apply(pd.to_numeric, errors="coerce").fillna(0)
    return frame


def read_obs(path: Path | None, *, h5ad_path: Path | None = None) -> pd.DataFrame | None:
    if path is not None:
        suffix = path.suffix.lower()
        if suffix in {".tsv", ".txt"}:
            return pd.read_csv(path, sep="\t", index_col=0)
        if suffix == ".csv":
            return pd.read_csv(path, index_col=0)
        if suffix in {".parquet", ".pq"}:
            return pd.read_parquet(path)
        raise ValueError(f"unsupported obs format: {path}")
    if h5ad_path is not None and h5ad_path.suffix.lower() == ".h5ad":
        obj = ad.read_h5ad(h5ad_path, backed="r")
        obs = obj.obs.copy()
        obs.index = obs.index.astype(str)
        obj.file.close()
        return obs
    return None


def write_outputs(
    *,
    result,
    output_dir: Path,
    sample_name: str,
    mode: str = "compact",
) -> dict[str, str]:
    if mode not in {"compact", "parquet", "tsv"}:
        raise ValueError("mode must be compact, parquet, or tsv")
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, str] = {}

    import json

    manifest_path = output_dir / f"{sample_name}_fastcnvpy_manifest.json"
    manifest_path.write_text(json.dumps(result.manifest, indent=2, sort_keys=True), encoding="utf-8")
    outputs["manifest"] = str(manifest_path)

    cell_meta_path = output_dir / f"{sample_name}_fastcnvpy_cell_metadata.tsv"
    result.cell_metadata.to_csv(cell_meta_path, sep="\t")
    outputs["cell_metadata"] = str(cell_meta_path)

    arm_path = output_dir / f"{sample_name}_fastcnvpy_arm_cnv.tsv"
    result.arm_cnv.to_csv(arm_path, sep="\t")
    outputs["arm_cnv"] = str(arm_path)

    windows_path = output_dir / f"{sample_name}_fastcnvpy_genomic_windows.tsv"
    result.window_metadata.to_csv(windows_path, sep="\t", index=False)
    outputs["genomic_windows"] = str(windows_path)

    if mode == "compact":
        return outputs

    if mode == "parquet":
        raw_path = output_dir / f"{sample_name}_rawGenomicScores.parquet"
        trim_path = output_dir / f"{sample_name}_genomicScores.parquet"
        result.raw_genomic_scores.to_parquet(raw_path)
        result.genomic_scores.to_parquet(trim_path)
    else:
        raw_path = output_dir / f"{sample_name}_rawGenomicScores.tsv"
        trim_path = output_dir / f"{sample_name}_genomicScores.tsv"
        result.raw_genomic_scores.to_csv(raw_path, sep="\t")
        result.genomic_scores.to_csv(trim_path, sep="\t")
    outputs["raw_genomic_scores"] = str(raw_path)
    outputs["genomic_scores"] = str(trim_path)
    return outputs


def write_pooled_outputs(
    *,
    result,
    output_dir: Path,
    sample_name: str,
    mode: str = "compact",
) -> dict[str, str]:
    if mode not in {"compact", "parquet", "tsv"}:
        raise ValueError("mode must be compact, parquet, or tsv")
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, str] = {}

    import json

    manifest_path = output_dir / f"{sample_name}_fastcnvpy_pooled_manifest.json"
    manifest_path.write_text(json.dumps(result.manifest, indent=2, sort_keys=True), encoding="utf-8")
    outputs["manifest"] = str(manifest_path)

    cell_meta_path = output_dir / f"{sample_name}_fastcnvpy_pooled_cell_metadata.tsv"
    result.cell_metadata.to_csv(cell_meta_path, sep="\t")
    outputs["cell_metadata"] = str(cell_meta_path)

    arm = result.arm_cnv.copy()
    if hasattr(arm.columns, "to_flat_index"):
        arm.columns = ["|".join(map(str, col)) if isinstance(col, tuple) else str(col) for col in arm.columns.to_flat_index()]
    arm_path = output_dir / f"{sample_name}_fastcnvpy_pooled_arm_cnv.tsv"
    arm.to_csv(arm_path, sep="\t")
    outputs["arm_cnv"] = str(arm_path)

    windows_path = output_dir / f"{sample_name}_fastcnvpy_pooled_genomic_windows.tsv"
    result.window_metadata.to_csv(windows_path, sep="\t", index=False)
    outputs["genomic_windows"] = str(windows_path)

    per_sample_dir = output_dir / "per_sample"
    per_sample_dir.mkdir(exist_ok=True)
    outputs["per_sample_dir"] = str(per_sample_dir)
    for sample, sample_result in result.per_sample.items():
        outputs[f"sample:{sample}"] = write_outputs(
            result=sample_result,
            output_dir=per_sample_dir / _safe_name(sample),
            sample_name=_safe_name(sample),
            mode=mode,
        )
    return outputs


def _safe_name(value: str) -> str:
    out = "".join(ch if ch.isalnum() or ch in {"-", "_", "."} else "_" for ch in str(value))
    return out or "sample"
