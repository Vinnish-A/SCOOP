from __future__ import annotations

from typing import Literal
import numpy as np
from scipy import sparse


def get_matrix(adata, layer: str | None = None):
    if layer is None:
        return adata.X
    if layer not in adata.layers:
        raise KeyError(f"Layer not found: {layer}")
    return adata.layers[layer]


def sum_axis(matrix, axis: Literal[0, 1]) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.sum(axis=axis)).ravel()
    return np.asarray(matrix.sum(axis=axis)).ravel()


def nnz_axis(matrix, axis: Literal[0, 1]) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.diff(matrix.tocsr().indptr) if axis == 1 else np.diff(matrix.tocsc().indptr)
    return np.count_nonzero(matrix, axis=axis)


def subset_gene_sum(adata, gene_mask, layer: str) -> np.ndarray:
    X = get_matrix(adata, layer)
    if gene_mask.sum() == 0:
        return np.zeros(adata.n_obs, dtype=float)
    return sum_axis(X[:, gene_mask], axis=1)
