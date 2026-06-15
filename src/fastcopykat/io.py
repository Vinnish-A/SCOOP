from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def read_counts(path: Path, *, layer: str | None = None, transpose: bool = False) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".h5ad":
        return read_h5ad_counts(path, layer=layer, transpose=transpose)
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
    return frame.apply(pd.to_numeric, errors="coerce").fillna(0)


def read_h5ad_counts(path: Path, *, layer: str | None = None, transpose: bool = False) -> pd.DataFrame:
    obj = ad.read_h5ad(path)
    matrix = obj.layers[layer] if layer else obj.X
    if sp.issparse(matrix):
        matrix = matrix.toarray()
    values = np.asarray(matrix)
    frame = pd.DataFrame(values, index=obj.obs_names.astype(str), columns=obj.var_names.astype(str))
    if not transpose:
        frame = frame.T
    return frame


def write_copykat_outputs(
    *,
    prediction: pd.DataFrame,
    cna: pd.DataFrame,
    cell_scores: pd.DataFrame | None = None,
    output_dir: Path,
    sample_name: str,
    manifest: dict,
    mode: str = "copykat-tsv",
) -> dict[str, str]:
    if mode not in {"compact", "parquet", "copykat-tsv"}:
        raise ValueError("mode must be one of: compact, parquet, copykat-tsv")
    output_dir.mkdir(parents=True, exist_ok=True)
    prediction_path = output_dir / f"{sample_name}_copykat_prediction.txt"
    manifest_path = output_dir / f"{sample_name}_fastcopykat_manifest.json"
    prediction.to_csv(prediction_path, sep="\t", index=False)

    import json

    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    outputs = {
        "prediction": str(prediction_path),
        "manifest": str(manifest_path),
    }
    if cell_scores is not None:
        scores_path = output_dir / f"{sample_name}_fastcopykat_cell_scores.tsv"
        cell_scores.to_csv(scores_path, sep="\t", index=False)
        outputs["cell_scores"] = str(scores_path)
    if mode == "compact":
        return outputs

    if mode == "parquet":
        cna_path = output_dir / f"{sample_name}_fastcopykat_CNA_final_results_bin_by_cell.parquet"
        cna.to_parquet(cna_path, index=False)
        outputs["cna"] = str(cna_path)
        return outputs

    cna_path = output_dir / f"{sample_name}_copykat_CNA_final_results_bin_by_cell.txt"
    cna.to_csv(cna_path, sep="\t", index=False)
    outputs["cna"] = str(cna_path)
    return outputs
