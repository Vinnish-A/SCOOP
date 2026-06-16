from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd

from .abundance_data import AbundanceTable, build_sample_celltype_counts_from_h5ad, read_sample_celltype_counts, write_abundance_inputs
from .abundance_design import merge_metadata
from .abundance_result import AbundanceResult
from .abundance_mil import run_mil_abundance


def load_abundance_table(
    input_h5ad: str | Path | None = None,
    counts: str | Path | None = None,
    sample_key: str = "sample_id",
    celltype_key: str = "cell_type_lvl3",
    min_cells_per_sample: int = 20,
    min_total_cells_per_celltype: int = 50,
) -> AbundanceTable:
    if input_h5ad is None and counts is None:
        raise ValueError("provide either input_h5ad or counts")
    if input_h5ad is not None:
        adata = ad.read_h5ad(input_h5ad)
        return build_sample_celltype_counts_from_h5ad(
            adata,
            sample_key=sample_key,
            celltype_key=celltype_key,
            min_cells_per_sample=min_cells_per_sample,
            min_total_cells_per_celltype=min_total_cells_per_celltype,
        )
    table = read_sample_celltype_counts(counts)
    table.sample_key = sample_key
    table.celltype_key = celltype_key
    return table


def _read_metadata(path: str | Path | None) -> pd.DataFrame | None:
    if path is None:
        return None
    meta = pd.read_csv(path, sep="\t")
    return meta


def _manifest(
    *,
    mode: str,
    table: AbundanceTable,
    sample_key: str,
    celltype_key: str,
    transform: str,
    n_events: int | None = None,
    n_classes: int | None = None,
    metrics: dict[str, Any] | None = None,
    model: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "schema_version": "fastde.abundance.v1",
        "engine": "fastde.abundance",
        "scsurvival_compat": True,
        "scsurvival_source": {
            "repository": "https://github.com/cliffren/scSurvival",
            "commit": "a76de0a00035e4c4d49df4a06001b392c8014105",
            "version": "1.3.0",
            "license": "GPL-3.0",
        },
        "mode": mode,
        "sample_key": sample_key,
        "celltype_key": celltype_key,
        "n_samples": table.n_samples,
        "n_celltypes": table.n_celltypes,
        "n_events": n_events,
        "n_classes": n_classes,
        "transform": transform,
        "model": model or {},
        "loss": {"task": mode},
        "metrics": metrics or {},
        "outputs": {},
        "diagnostics": {},
    }


def run_abundance(
    *,
    mode: str,
    output_dir: str | Path,
    input_h5ad: str | Path | None = None,
    counts: str | Path | None = None,
    metadata: str | Path | pd.DataFrame | None = None,
    sample_key: str = "sample_id",
    celltype_key: str = "cell_type_lvl3",
    label_col: str | None = None,
    positive_label: str | None = None,
    negative_label: str | None = None,
    reference_level: str | None = None,
    time_col: str | None = None,
    event_col: str | None = None,
    value_col: str | None = None,
    covariates: str | list[str] | None = None,
    min_cells_per_sample: int = 20,
    min_total_cells_per_celltype: int = 50,
    **kwargs: Any,
) -> AbundanceResult:
    mode = mode.lower()
    if mode == "condition":
        mode = "binary"
    table = load_abundance_table(
        input_h5ad=input_h5ad,
        counts=counts,
        sample_key=sample_key,
        celltype_key=celltype_key,
        min_cells_per_sample=min_cells_per_sample,
        min_total_cells_per_celltype=min_total_cells_per_celltype,
    )
    if isinstance(metadata, pd.DataFrame):
        meta = metadata
    else:
        meta = _read_metadata(metadata)
    table = merge_metadata(table, meta, sample_key=sample_key)
    outdir = Path(output_dir)
    write_abundance_inputs(table, outdir)

    return run_mil_abundance(
        mode=mode,
        table=table,
        output_dir=outdir,
        input_h5ad=input_h5ad,
        sample_key=sample_key,
        celltype_key=celltype_key,
        label_col=label_col,
        positive_label=positive_label,
        negative_label=negative_label,
        reference_level=reference_level,
        time_col=time_col,
        event_col=event_col,
        value_col=value_col,
        covariates=covariates,
        max_instances_per_sample=int(kwargs.get("max_instances_per_sample", 2000)),
        hidden_dim=int(kwargs.get("hidden_dim", 64)),
        dropout=float(kwargs.get("dropout", 0.1)),
        learning_rate=float(kwargs.get("learning_rate", 1e-3)),
        weight_decay=float(kwargs.get("weight_decay", 1e-4)),
        max_epochs=int(kwargs.get("max_epochs", 500)),
        random_seed=int(kwargs.get("random_seed", 0)),
        survival_loss=str(kwargs.get("survival_loss", "cox")),
        manifest_factory=_manifest,
    )
