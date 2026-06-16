from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class AbundanceTable:
    counts: pd.DataFrame
    proportions: pd.DataFrame
    metadata: pd.DataFrame
    sample_key: str = "sample_id"
    celltype_key: str = "cell_type"

    @property
    def n_samples(self) -> int:
        return int(self.counts.shape[0])

    @property
    def n_celltypes(self) -> int:
        return int(self.counts.shape[1])


def _proportions(counts: pd.DataFrame) -> pd.DataFrame:
    totals = counts.sum(axis=1).replace(0, np.nan)
    return counts.div(totals, axis=0).fillna(0.0)


def build_sample_celltype_counts_from_h5ad(
    adata,
    sample_key: str,
    celltype_key: str,
    min_cells_per_sample: int = 20,
    min_total_cells_per_celltype: int = 50,
) -> AbundanceTable:
    if sample_key not in adata.obs:
        raise KeyError(f"sample_key {sample_key!r} is not present in adata.obs")
    if celltype_key not in adata.obs:
        raise KeyError(f"celltype_key {celltype_key!r} is not present in adata.obs")

    obs = adata.obs[[sample_key, celltype_key]].copy()
    obs = obs.dropna()
    counts = pd.crosstab(obs[sample_key].astype(str), obs[celltype_key].astype(str))
    counts = counts.sort_index(axis=0).sort_index(axis=1).astype(int)

    sample_totals = counts.sum(axis=1)
    counts = counts.loc[sample_totals >= min_cells_per_sample]
    celltype_totals = counts.sum(axis=0)
    counts = counts.loc[:, celltype_totals >= min_total_cells_per_celltype]
    if counts.empty:
        raise ValueError("no samples or cell types remain after abundance filtering")

    metadata = pd.DataFrame(index=counts.index)
    metadata.index.name = sample_key
    return AbundanceTable(
        counts=counts,
        proportions=_proportions(counts),
        metadata=metadata,
        sample_key=sample_key,
        celltype_key=celltype_key,
    )


def read_sample_celltype_counts(path: str | Path) -> AbundanceTable:
    counts = pd.read_csv(path, sep="\t", index_col=0)
    counts.index = counts.index.astype(str)
    counts.columns = counts.columns.astype(str)
    counts = counts.fillna(0).astype(float)
    metadata = pd.DataFrame(index=counts.index)
    metadata.index.name = "sample_id"
    return AbundanceTable(counts=counts, proportions=_proportions(counts), metadata=metadata)


def write_abundance_inputs(table: AbundanceTable, output_dir: str | Path) -> dict[str, str]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    counts_path = out / "sample_by_celltype_counts.tsv"
    prop_path = out / "sample_by_celltype_proportions.tsv"
    meta_path = out / "abundance_metadata_used.tsv"
    table.counts.to_csv(counts_path, sep="\t")
    table.proportions.to_csv(prop_path, sep="\t")
    table.metadata.to_csv(meta_path, sep="\t")
    return {
        "counts": str(counts_path),
        "proportions": str(prop_path),
        "metadata": str(meta_path),
    }
