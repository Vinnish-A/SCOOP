"""Compatibility utilities for using harmonypy 2.x with cNMF-style MOE correction.

harmonypy 2.x exposes the corrected PCA embedding through a C++ backend, but it
does not expose the ``Phi_moe`` design matrix or the final ridge lambda vector
that cNMF 1.7.x expects when applying Harmony's MOE correction to the normalized
gene matrix. This module reconstructs those quantities for the fixed-lambda
case and applies the same ridge correction formula used by cNMF.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class Harmony2MOEResult:
    """Result returned by ``harmony2_moe_correct``."""

    x_corr: np.ndarray
    x_pca_harmony: np.ndarray
    phi_moe: np.ndarray
    lamb: np.ndarray
    r: np.ndarray
    k: int


def normalize_vars_use(vars_use: str | Iterable[str]) -> list[str]:
    if isinstance(vars_use, str):
        return [vars_use]
    return list(vars_use)


def build_phi_moe(obs, vars_use: str | Iterable[str]) -> tuple[np.ndarray, list[int]]:
    """Rebuild Harmony's one-hot MOE design matrix with intercept.

    The returned matrix has shape ``(1 + n_batch_levels, n_cells)``. Category
    ordering follows ``np.unique``, matching harmonypy 2.0's compact
    ``batch_of_cell`` construction.
    """

    vars_use = normalize_vars_use(vars_use)
    n_cells = obs.shape[0]
    blocks: list[np.ndarray] = []
    phi_n: list[int] = []
    for var in vars_use:
        values = np.asarray(obs[var])
        _, codes = np.unique(values, return_inverse=True)
        n_levels = int(codes.max() + 1)
        phi_n.append(n_levels)
        block = np.zeros((n_levels, n_cells), dtype=np.float64)
        block[codes, np.arange(n_cells)] = 1.0
        blocks.append(block)
    phi = np.vstack(blocks)
    intercept = np.ones((1, n_cells), dtype=np.float64)
    return np.vstack([intercept, phi]), phi_n


def build_fixed_lamb(lamb: float | int | Iterable[float], phi_n: list[int]) -> np.ndarray:
    """Build the cNMF/Harmony fixed-lambda vector.

    cNMF 1.7.x passes this vector directly into ``gram + lamb``. NumPy then
    broadcasts the vector across Gram-matrix columns. FastCNMF keeps that
    behavior for legacy-compatible outputs.
    """

    if isinstance(lamb, (float, int)):
        values = np.repeat(float(lamb), sum(phi_n))
    else:
        arr = np.asarray(list(lamb), dtype=np.float64)
        if arr.size == len(phi_n):
            values = np.repeat(arr, phi_n)
        elif arr.size == sum(phi_n):
            values = arr
        elif arr.size == 1:
            values = np.repeat(float(arr[0]), sum(phi_n))
        else:
            raise ValueError(
                f"lamb length {arr.size} is incompatible with {len(phi_n)} "
                f"batch covariates and {sum(phi_n)} total levels"
            )
    return np.insert(values.astype(np.float64), 0, 0.0)


def normalize_lamb_addend(lamb: np.ndarray | Iterable[float]) -> np.ndarray:
    """Normalize Harmony lambda to the addend layout used by cNMF."""

    arr = np.asarray(lamb, dtype=np.float64)
    if arr.ndim == 1:
        return arr
    if arr.ndim == 2 and arr.shape[0] == arr.shape[1]:
        return arr
    raise ValueError(f"lamb must be a vector or square matrix, got shape {arr.shape}")


def moe_correct_ridge_fast(
    x_cells_by_genes: np.ndarray,
    r_cells_by_clusters: np.ndarray,
    phi_moe_features_by_cells: np.ndarray,
    lamb_vector: np.ndarray,
) -> np.ndarray:
    """Apply cNMF's Harmony MOE ridge correction to a cell x gene matrix."""

    z_orig = np.asarray(x_cells_by_genes, dtype=np.float64).T
    r = np.asarray(r_cells_by_clusters, dtype=np.float64).T
    phi_moe = np.asarray(phi_moe_features_by_cells, dtype=np.float64)
    lamb_addend = normalize_lamb_addend(lamb_vector)

    z_corr = z_orig.copy()
    for cluster_i in range(r.shape[0]):
        phi_rk = phi_moe * r[cluster_i, :]
        gram = phi_rk @ phi_moe.T + lamb_addend
        weights = np.linalg.solve(gram, phi_rk @ z_orig.T)
        weights[0, :] = 0.0
        z_corr -= weights.T @ phi_rk
    out = z_corr.T
    out[out < 0] = 0.0
    return out


def moe_correct_ridge_batched(
    x_cells_by_genes: np.ndarray,
    r_cells_by_clusters: np.ndarray,
    phi_moe_features_by_cells: np.ndarray,
    lamb_vector: np.ndarray,
) -> np.ndarray:
    """Batched MOE ridge correction optimized for one-hot Harmony covariates.

    This is algebraically equivalent to ``moe_correct_ridge_fast`` but avoids
    doing two small-by-huge matrix multiplications for each Harmony cluster.
    Instead, it batches the cluster dimension into a small number of large BLAS
    calls, one for each MOE design row.
    """

    x = np.asarray(x_cells_by_genes, dtype=np.float64)
    r = np.asarray(r_cells_by_clusters, dtype=np.float64)
    phi = np.asarray(phi_moe_features_by_cells, dtype=np.float64)
    lamb = normalize_lamb_addend(lamb_vector)

    n_cells, n_genes = x.shape
    n_features, phi_cells = phi.shape
    if phi_cells != n_cells:
        raise ValueError(f"phi_moe has {phi_cells} cells but X has {n_cells}")
    if r.shape[0] != n_cells:
        raise ValueError(f"R has {r.shape[0]} cells but X has {n_cells}")

    n_clusters = r.shape[1]
    # b_all[k, p, g] = sum_n R[n,k] * Phi[p,n] * X[n,g]
    b_all = np.empty((n_clusters, n_features, n_genes), dtype=np.float64)
    for p in range(n_features):
        b_all[:, p, :] = (r * phi[p, :, None]).T @ x

    # a_all[k, p, q] = sum_n R[n,k] * Phi[p,n] * Phi[q,n] + lamb[p,q]
    a_all = np.empty((n_clusters, n_features, n_features), dtype=np.float64)
    for p in range(n_features):
        for q in range(n_features):
            a_all[:, p, q] = r.T @ (phi[p, :] * phi[q, :])
    if lamb.ndim == 1:
        a_all += lamb[None, None, :]
    else:
        a_all += lamb[None, :, :]

    coef = np.linalg.solve(a_all, b_all)
    coef[:, 0, :] = 0.0

    correction = np.zeros_like(x)
    for p in range(1, n_features):
        correction += (r * phi[p, :, None]) @ coef[:, p, :]

    out = x - correction
    out[out < 0] = 0.0
    return out


def harmony2_moe_correct(
    x_cells_by_genes: np.ndarray,
    pca_cells_by_components: np.ndarray,
    obs,
    vars_use: str | Iterable[str],
    *,
    lamb: float | int | Iterable[float] = 1.0,
    theta: float | int | Iterable[float] = 1.0,
    max_iter_harmony: int = 20,
    ncores: int = 1,
    verbose: bool = False,
    random_state: int = 0,
    correction_method: str = "batched",
):
    """Run harmonypy 2.x then apply cNMF-compatible MOE gene correction.

    This function intentionally supports only fixed-lambda mode. Dynamic lambda
    compatibility requires harmonypy 2.x to expose the final lambda vector from
    its C++ backend.
    """

    import harmonypy

    phi_moe, phi_n = build_phi_moe(obs, vars_use)
    lamb_vector = build_fixed_lamb(lamb, phi_n)
    harmony_res = harmonypy.run_harmony(
        pca_cells_by_components,
        obs,
        vars_use,
        lamb=lamb,
        theta=theta,
        max_iter_harmony=max_iter_harmony,
        ncores=ncores,
        verbose=verbose,
        random_state=random_state,
    )
    r = np.asarray(harmony_res.R, dtype=np.float64)
    if correction_method == "batched":
        x_corr = moe_correct_ridge_batched(x_cells_by_genes, r, phi_moe, lamb_vector)
    elif correction_method == "cluster_loop":
        x_corr = moe_correct_ridge_fast(x_cells_by_genes, r, phi_moe, lamb_vector)
    else:
        raise ValueError(f"Unknown correction_method: {correction_method}")
    return Harmony2MOEResult(
        x_corr=x_corr,
        x_pca_harmony=np.asarray(harmony_res.Z_corr, dtype=np.float64),
        phi_moe=phi_moe,
        lamb=lamb_vector,
        r=r,
        k=int(harmony_res.K),
    )
