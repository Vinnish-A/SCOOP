from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from fastde.abundance_data import build_sample_celltype_counts_from_h5ad
from fastde.abundance_design import abundance_transform


def test_build_sample_celltype_counts_from_h5ad() -> None:
    obs = pd.DataFrame(
        {
            "sample_id": ["S1", "S1", "S1", "S2", "S2", "S3"],
            "cell_type_lvl3": ["T", "T", "B", "T", "M", "M"],
        }
    )
    adata = AnnData(np.ones((len(obs), 2)), obs=obs)
    table = build_sample_celltype_counts_from_h5ad(adata, "sample_id", "cell_type_lvl3", min_cells_per_sample=2, min_total_cells_per_celltype=1)
    assert table.counts.loc["S1", "T"] == 2
    assert table.counts.loc["S2", "B"] == 0
    assert "S3" not in table.counts.index


def test_clr_transform() -> None:
    counts = pd.DataFrame([[10, 0, 5], [2, 8, 0]], index=["S1", "S2"], columns=["A", "B", "C"])
    z = abundance_transform(counts, transform="clr", pseudocount=0.5)
    assert np.allclose(z.mean(axis=1), 0.0)
    assert np.isfinite(z.to_numpy()).all()
