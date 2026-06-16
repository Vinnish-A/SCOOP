from __future__ import annotations

from pathlib import Path
import shlex
import subprocess
from typing import Any, Mapping

import pandas as pd


def validate_lr_resource(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Ligand-receptor resource not found: {path}")
    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
    else:
        df = pd.read_csv(path, sep="\t" if path.suffix in {".tsv", ".txt"} else ",")
    required = {"ligand", "receptor"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"LR resource missing required columns: {sorted(missing)}")
    for col in ["ligand_subunits", "receptor_subunits"]:
        if col not in df.columns:
            df[col] = df[col.replace("_subunits", "")].astype(str)
    df["ligand_n_subunits"] = df["ligand_subunits"].astype(str).str.split(r"[|;,]").map(len)
    df["receptor_n_subunits"] = df["receptor_subunits"].astype(str).str.split(r"[|;,]").map(len)
    df["complex_sensitive"] = (df["ligand_n_subunits"] > 1) | (df["receptor_n_subunits"] > 1)
    return df


def choose_ccc_groupby(adata, primary: str, fallback: str, min_cells_per_group: int = 20) -> str:
    if primary in adata.obs:
        vc = adata.obs[primary].value_counts()
        if len(vc) and vc.min() >= min_cells_per_group:
            return primary
    if fallback in adata.obs:
        return fallback
    raise KeyError(f"No usable CCC groupby key: {primary} or {fallback}")


def build_fastccc_command(command_template: str, **kwargs) -> list[str]:
    return shlex.split(command_template.format(**kwargs))


def run_fastccc_command(command: list[str], cwd: str | Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(command, cwd=cwd, check=True, text=True, capture_output=True)


def run_cellphonedb_validation_omicverse(
    adata,
    *,
    cpdb_file_path: str,
    celltype_key: str,
    output_dir: str | Path,
    iterations: int = 1000,
    threshold: float = 0.1,
    pvalue: float = 0.05,
    threads: int = 10,
) -> dict[str, Any]:
    """Run OmicVerse CellPhoneDB wrapper and immediately externalise results."""
    from omicverse_transfer.external import require_omicverse

    ov = require_omicverse()
    results_key = "_cpdb_results_tmp"
    comm_key = "_cpdb_comm_tmp"
    raw, comm = ov.single.run_cellphonedb_v5(
        adata,
        cpdb_file_path=cpdb_file_path,
        celltype_key=celltype_key,
        iterations=iterations,
        threshold=threshold,
        pvalue=pvalue,
        threads=threads,
        output_dir=str(output_dir),
        cleanup_temp=True,
        debug=False,
        results_key=results_key,
        comm_key=comm_key,
    )
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    exported = {}
    if isinstance(raw, dict):
        for name, obj in raw.items():
            if isinstance(obj, pd.DataFrame):
                path = out_dir / f"cpdb_{name}.parquet"
                obj.to_parquet(path, index=False)
                exported[name] = str(path)
    # Keep H5AD lean.
    for key in [results_key, comm_key]:
        if key in adata.uns:
            del adata.uns[key]
    return {"exported": exported, "comm_n_obs": getattr(comm, "n_obs", None)}


def run_liana_validation_omicverse(adata, *, groupby: str, method: str = "rank_aggregate") -> pd.DataFrame:
    from omicverse_transfer.external import require_omicverse

    ov = require_omicverse()
    res = ov.single.run_liana(adata, groupby=groupby, method=method, inplace=False)
    if not isinstance(res, pd.DataFrame):
        raise TypeError("OmicVerse run_liana did not return a DataFrame with inplace=False")
    return res
