from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

from fastde.deseq2 import estimate_size_factors, run_deseq2_wald
from fastde.markers import run_cosg_markers, run_cosg_markers_dense, run_cosg_markers_sparse


def test_deseq2_wald_detects_up_direction() -> None:
    rng = np.random.default_rng(1)
    counts = pd.DataFrame(
        rng.poisson(80, size=(8, 20)),
        index=[f"s{i}" for i in range(8)],
        columns=[f"G{i}" for i in range(20)],
    )
    counts.loc[[f"s{i}" for i in range(4, 8)], "G0"] *= 4
    meta = pd.DataFrame({"condition": ["ctrl"] * 4 + ["test"] * 4}, index=counts.index)

    result = run_deseq2_wald(counts, meta, condition_col="condition", ctrl_group="ctrl", test_group="test")
    table = result.table.set_index("gene")

    assert table.loc["G0", "log2FoldChange"] > 1.0
    assert table.loc["G0", "pvalue"] < 0.05
    assert {"dispGeneEst", "dispFit", "dispMAP", "dispOutlier"}.issubset(result.table.columns)
    assert result.manifest["method"] == "fastde_deseq2_wald"
    assert result.manifest["dispersion_trend"]["type"] in {"parametric", "mean", "mean_fallback"}


def test_size_factors_are_centered_geometrically() -> None:
    counts = np.array([[10, 20, 30], [20, 40, 60], [5, 10, 15]], dtype=float)
    sf = estimate_size_factors(counts)
    assert np.isclose(np.exp(np.mean(np.log(sf))), 1.0)
    assert sf[1] > sf[0] > sf[2]


def test_cosg_markers_prioritize_group_specific_gene() -> None:
    matrix = np.array(
        [
            [5, 0, 1],
            [6, 0, 1],
            [0, 5, 1],
            [0, 6, 1],
        ],
        dtype=float,
    )
    labels = ["A", "A", "B", "B"]
    result = run_cosg_markers(matrix, labels, ["GA", "GB", "HOUSE"], top_n=2)
    top = result.table[result.table["rank"] == 1].set_index("group")

    assert top.loc["A", "gene"] == "GA"
    assert top.loc["B", "gene"] == "GB"


def test_sparse_cosg_matches_dense_reference() -> None:
    matrix = np.array(
        [
            [5, 0, 1, 0],
            [6, 0, 1, 0],
            [0, 5, 1, 2],
            [0, 6, 1, 2],
            [1, 0, 0, 0],
            [0, 1, 0, 0],
        ],
        dtype=float,
    )
    labels = ["A", "A", "B", "B", "A", "B"]
    genes = ["GA", "GB", "HOUSE", "B2"]

    dense = run_cosg_markers_dense(matrix, labels, genes, top_n=3, min_pct=0.0).table
    sparse_result = run_cosg_markers_sparse(sparse.csr_matrix(matrix), labels, genes, top_n=3, min_pct=0.0).table

    dense_key = dense.set_index(["group", "rank"])
    sparse_key = sparse_result.set_index(["group", "rank"])
    assert dense_key["gene"].to_dict() == sparse_key["gene"].to_dict()
    assert np.allclose(dense_key["score"], sparse_key["score"])
