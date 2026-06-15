from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

import numpy as np
import pandas as pd
from scipy import sparse, stats


@dataclass(frozen=True)
class MarkerResult:
    table: pd.DataFrame
    manifest: dict


def run_cosg_markers(
    matrix,
    labels,
    gene_names,
    *,
    groups: list[str] | None = None,
    top_n: int = 100,
    min_pct: float = 0.05,
) -> MarkerResult:
    if _is_sparse_like(matrix):
        return run_cosg_markers_sparse(matrix, labels, gene_names, groups=groups, top_n=top_n, min_pct=min_pct)
    return run_cosg_markers_dense(matrix, labels, gene_names, groups=groups, top_n=top_n, min_pct=min_pct)


def run_cosg_markers_dense(
    matrix,
    labels,
    gene_names,
    *,
    groups: list[str] | None = None,
    top_n: int = 100,
    min_pct: float = 0.05,
) -> MarkerResult:
    start = perf_counter()
    x = _as_dense_float(matrix)
    labels = pd.Series(labels).astype(str).to_numpy()
    genes = pd.Index(gene_names).astype(str)
    groups = groups or sorted(pd.unique(labels))
    x_norm = np.sqrt((x * x).sum(axis=0))
    x_norm = np.where(x_norm > 0, x_norm, np.inf)
    rows = []
    expressed = x > 0
    for group in groups:
        mask = labels == str(group)
        if mask.sum() == 0:
            continue
        y_norm = np.sqrt(float(mask.sum()))
        score = x[mask, :].sum(axis=0) / (x_norm * y_norm)
        pct_in = expressed[mask, :].mean(axis=0)
        pct_out = expressed[~mask, :].mean(axis=0) if (~mask).sum() else np.zeros(x.shape[1])
        mean_in = x[mask, :].mean(axis=0)
        mean_out = x[~mask, :].mean(axis=0) if (~mask).sum() else np.zeros(x.shape[1])
        keep = pct_in >= min_pct
        order = np.lexsort((genes.to_numpy(), -score))
        selected = [idx for idx in order if keep[idx]][:top_n]
        for rank, idx in enumerate(selected, start=1):
            rows.append(
                {
                    "group": group,
                    "rank": rank,
                    "gene": genes[idx],
                    "score": float(score[idx]),
                    "pct_in": float(pct_in[idx]),
                    "pct_out": float(pct_out[idx]),
                    "mean_in": float(mean_in[idx]),
                    "mean_out": float(mean_out[idx]),
                    "mean_diff": float(mean_in[idx] - mean_out[idx]),
                }
            )
    table = pd.DataFrame(rows)
    manifest = {
        "schema_version": "fastde.markers.cosg_like.v1",
        "method": "cosg_like",
        "n_cells": int(x.shape[0]),
        "n_genes": int(x.shape[1]),
        "groups": list(map(str, groups)),
        "top_n": int(top_n),
        "min_pct": float(min_pct),
        "backend": "dense_reference",
        "seconds": round(perf_counter() - start, 6),
    }
    return MarkerResult(table=table, manifest=manifest)


def run_cosg_markers_sparse(
    matrix,
    labels,
    gene_names,
    *,
    groups: list[str] | None = None,
    top_n: int = 100,
    min_pct: float = 0.05,
) -> MarkerResult:
    """COSG-like marker scoring without densifying the cell-by-gene matrix.

    The core score is the same as :func:`run_cosg_markers_dense`:
    cosine similarity between each gene vector and a binary group-membership
    vector. The sparse implementation computes all group sums as ``G @ X``,
    where ``G`` is a group-by-cell sparse indicator matrix, and only densifies
    the small group-by-gene summaries.
    """

    start = perf_counter()
    x = _as_csc_float(matrix)
    labels = pd.Series(labels).astype(str).to_numpy()
    if x.shape[0] != labels.size:
        raise ValueError("matrix row count must match labels length")
    genes = pd.Index(gene_names).astype(str)
    if x.shape[1] != len(genes):
        raise ValueError("matrix column count must match gene_names length")

    groups = groups or sorted(pd.unique(labels))
    groups = list(map(str, groups))
    group_to_code = {group: i for i, group in enumerate(groups)}
    codes = np.array([group_to_code.get(str(label), -1) for label in labels], dtype=np.int64)
    valid = codes >= 0
    if not np.any(valid):
        raise ValueError("none of the requested groups are present in labels")

    group_codes = codes[valid]
    group_n = np.bincount(group_codes, minlength=len(groups)).astype(np.float64)
    group_sum, group_nnz, total_sum, total_nnz, x_norm = _group_summaries_by_csc_column(x, codes, len(groups))
    x_norm = np.where(x_norm > 0, x_norm, np.inf)

    rows = []
    gene_values = genes.to_numpy()
    n_cells = float(x.shape[0])
    for g, group in enumerate(groups):
        n_in = float(group_n[g])
        if n_in == 0:
            continue
        n_out = n_cells - n_in
        score = group_sum[g] / (x_norm * np.sqrt(n_in))
        pct_in = group_nnz[g] / n_in
        pct_out = (total_nnz - group_nnz[g]) / n_out if n_out > 0 else np.zeros(x.shape[1], dtype=float)
        mean_in = group_sum[g] / n_in
        mean_out = (total_sum - group_sum[g]) / n_out if n_out > 0 else np.zeros(x.shape[1], dtype=float)
        keep = pct_in >= min_pct
        order = np.lexsort((gene_values, -score))
        selected = order[keep[order]][:top_n]
        for rank, idx in enumerate(selected, start=1):
            rows.append(
                {
                    "group": group,
                    "rank": rank,
                    "gene": genes[idx],
                    "score": float(score[idx]),
                    "pct_in": float(pct_in[idx]),
                    "pct_out": float(pct_out[idx]),
                    "mean_in": float(mean_in[idx]),
                    "mean_out": float(mean_out[idx]),
                    "mean_diff": float(mean_in[idx] - mean_out[idx]),
                }
            )
    table = pd.DataFrame(rows)
    manifest = {
        "schema_version": "fastde.markers.cosg_like.v1",
        "method": "cosg_like",
        "n_cells": int(x.shape[0]),
        "n_genes": int(x.shape[1]),
        "n_nonzero": int(x.nnz),
        "groups": groups,
        "top_n": int(top_n),
        "min_pct": float(min_pct),
        "backend": "sparse_csc_column_scan",
        "summary_matrix_shape": [int(len(groups)), int(x.shape[1])],
        "seconds": round(perf_counter() - start, 6),
    }
    return MarkerResult(table=table, manifest=manifest)


def run_wilcoxon_markers(
    matrix,
    labels,
    gene_names,
    *,
    groups: list[str] | None = None,
    top_n: int = 100,
    min_pct: float = 0.05,
) -> MarkerResult:
    start = perf_counter()
    x = _as_dense_float(matrix)
    labels = pd.Series(labels).astype(str).to_numpy()
    genes = pd.Index(gene_names).astype(str)
    groups = groups or sorted(pd.unique(labels))
    rows = []
    expressed = x > 0
    for group in groups:
        mask = labels == str(group)
        pct_in = expressed[mask, :].mean(axis=0)
        pct_out = expressed[~mask, :].mean(axis=0) if (~mask).sum() else np.zeros(x.shape[1])
        mean_in = x[mask, :].mean(axis=0)
        mean_out = x[~mask, :].mean(axis=0) if (~mask).sum() else np.zeros(x.shape[1])
        pvals = np.ones(x.shape[1], dtype=float)
        stats_val = np.zeros(x.shape[1], dtype=float)
        for j in range(x.shape[1]):
            if pct_in[j] < min_pct:
                continue
            res = stats.ranksums(x[mask, j], x[~mask, j])
            stats_val[j] = res.statistic
            pvals[j] = res.pvalue
        padj = benjamini_hochberg(pvals)
        order = np.lexsort((genes.to_numpy(), pvals, -stats_val))
        selected = [idx for idx in order if pct_in[idx] >= min_pct][:top_n]
        for rank, idx in enumerate(selected, start=1):
            rows.append(
                {
                    "group": group,
                    "rank": rank,
                    "gene": genes[idx],
                    "score": float(stats_val[idx]),
                    "pvalue": float(pvals[idx]),
                    "padj": float(padj[idx]),
                    "pct_in": float(pct_in[idx]),
                    "pct_out": float(pct_out[idx]),
                    "mean_in": float(mean_in[idx]),
                    "mean_out": float(mean_out[idx]),
                    "mean_diff": float(mean_in[idx] - mean_out[idx]),
                }
            )
    table = pd.DataFrame(rows)
    manifest = {
        "schema_version": "fastde.markers.wilcoxon.v1",
        "method": "wilcoxon",
        "n_cells": int(x.shape[0]),
        "n_genes": int(x.shape[1]),
        "groups": list(map(str, groups)),
        "top_n": int(top_n),
        "min_pct": float(min_pct),
        "seconds": round(perf_counter() - start, 6),
    }
    return MarkerResult(table=table, manifest=manifest)


def _as_dense_float(matrix) -> np.ndarray:
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=np.float64)


def _is_sparse_like(matrix) -> bool:
    return sparse.issparse(matrix) or hasattr(matrix, "to_memory")


def _as_csc_float(matrix) -> sparse.csc_matrix:
    if hasattr(matrix, "to_memory") and not sparse.issparse(matrix):
        matrix = matrix.to_memory()
    if not sparse.issparse(matrix):
        raise TypeError("matrix is not sparse")
    if not sparse.isspmatrix_csc(matrix):
        matrix = matrix.tocsc()
    return matrix


def _group_summaries_by_csc_column(x: sparse.csc_matrix, codes: np.ndarray, n_groups: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n_genes = x.shape[1]
    group_sum = np.zeros((n_groups, n_genes), dtype=np.float64)
    group_nnz = np.zeros((n_groups, n_genes), dtype=np.float64)
    total_sum = np.zeros(n_genes, dtype=np.float64)
    total_nnz = np.zeros(n_genes, dtype=np.float64)
    sq_sum = np.zeros(n_genes, dtype=np.float64)
    indptr = x.indptr
    indices = x.indices
    data = x.data
    for j in range(n_genes):
        start, end = indptr[j], indptr[j + 1]
        if start == end:
            continue
        rows = indices[start:end]
        values = np.asarray(data[start:end], dtype=np.float64)
        total_sum[j] = values.sum()
        total_nnz[j] = float(values.size)
        sq_sum[j] = float(np.dot(values, values))
        group_codes = codes[rows]
        valid = group_codes >= 0
        if not np.any(valid):
            continue
        valid_codes = group_codes[valid]
        group_sum[:, j] = np.bincount(valid_codes, weights=values[valid], minlength=n_groups)
        group_nnz[:, j] = np.bincount(valid_codes, minlength=n_groups)
    return group_sum, group_nnz, total_sum, total_nnz, np.sqrt(sq_sum)


def benjamini_hochberg(pvalue: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalue, dtype=float)
    n = p.size
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out = np.empty_like(adjusted)
    out[order] = adjusted
    return out
