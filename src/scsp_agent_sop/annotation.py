from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm


def run_markers_scanpy(adata, groupby: str = "cluster_identity", layer: str = "log1p_norm", method: str = "wilcoxon",
                       n_genes: int | None = None, measurement=None, min_measured: int = 3) -> pd.DataFrame:
    """One-vs-rest marker tests; with ``measurement``, partially measured genes use measuring cells only.

    Genes measured in every sample of ``adata`` go through scanpy unchanged. The
    other genes are grouped by measured-sample pattern and tested with an exact
    tie-corrected two-sided Mann-Whitney test (scanpy's wilcoxon statistic) on
    the cells whose sample measured them. A group with fewer than
    ``min_measured`` measuring cells on either side gets no row for that gene
    (unassessed). Wilcoxon output adds ``U``, ``n_measured_group``,
    ``n_measured_rest``, ``log10_pval`` and ``log10_pval_adj`` (BH in log10
    space), which stay finite where ``pvals`` underflows to 0.
    """
    import scanpy as sc
    if layer != "X" and layer not in adata.layers:
        raise ValueError(f"Requested differential-expression layer is missing: {layer}")
    full = None
    if measurement is not None:
        measurement.check(adata)
        full = measurement.measured_in_all(adata.obs[measurement.sample_key])
    masked = full is not None and not full.all()
    if masked and method != "wilcoxon":
        raise ValueError("Partially measured genes are tested with a masked Mann-Whitney; use method='wilcoxon'")
    sc.tl.rank_genes_groups(adata, groupby=groupby, layer=None if layer == "X" else layer,
                           use_raw=False, method=method, n_genes=None if masked else n_genes, key_added="_markers_tmp",
                           tie_correct=True, pts=True, mask_var=full if masked else None)
    df = sc.get.rank_genes_groups_df(adata, group=None, key="_markers_tmp")
    if "_markers_tmp" in adata.uns:
        del adata.uns["_markers_tmp"]
    if method != "wilcoxon":
        return df

    matrix = adata.X if layer == "X" else adata.layers[layer]
    groups = pd.Categorical(adata.obs[groupby])
    categories = np.asarray(groups.categories).astype(str)
    # cells without a group sit in the last bin: scanpy counts them in "rest"
    bins = np.where(groups.codes < 0, len(categories), groups.codes)
    full_genes = np.flatnonzero(full) if masked else np.arange(adata.n_vars)
    stats = [_mwu_rows(matrix, None, full_genes, bins, categories, adata.var_names, 0)]
    if masked:
        codes = measurement.sample_codes(adata.obs[measurement.sample_key])
        present = np.unique(codes)
        partial = np.flatnonzero(~full)
        patterns, pattern_of_gene = np.unique(measurement.measured[np.ix_(present, partial)].T, axis=0, return_inverse=True)
        for p, pattern in enumerate(patterns):
            if pattern.any():
                cells = np.flatnonzero(np.isin(codes, present[pattern]))
                stats.append(_mwu_rows(matrix, cells, partial[pattern_of_gene.ravel() == p], bins[cells],
                                       categories, adata.var_names, min_measured))
    stats = pd.concat(stats, ignore_index=True)
    stats["log10_pval_adj"] = stats.groupby("group").log10_pval.transform(lambda s: _log10_bh(s.to_numpy()))
    extra = ["group", "names", "n_measured_group", "n_measured_rest", "U", "log10_pval", "log10_pval_adj"]
    df = df.merge(stats[extra], on=["group", "names"], how="left")
    if not masked:
        return df
    tested = stats[~stats.names.isin(df.names)].copy()
    tested["pvals"] = 2 * norm.sf(np.abs(tested.scores.to_numpy()))
    df = pd.concat([df, tested], ignore_index=True)
    df["pvals_adj"] = 10.0 ** df.log10_pval_adj  # one BH family per group over both routes
    df["_group_order"] = pd.Categorical(df.group, categories=categories).codes
    df = df.sort_values(["_group_order", "scores"], ascending=[True, False], kind="mergesort").drop(columns="_group_order")
    if n_genes is not None:
        df = df.groupby("group", sort=False).head(n_genes)
    return df.reset_index(drop=True)


def _mwu_rows(matrix, cells, genes, bins, categories, var_names, min_measured: int, chunk: int = 1024) -> pd.DataFrame:
    """Exact tie-corrected one-vs-rest Mann-Whitney over ``cells`` (all when None) for ``genes``.

    Stored entries are ranked with one tied block for the implicit zeros, so a
    sparse column is never densified.
    """
    import scipy.sparse as sp

    k = len(categories)
    n_bins = k + 1
    n = len(bins)
    bin_counts = np.bincount(bins, minlength=n_bins)
    n_group = bin_counts[:k].astype(np.float64)
    shape = (n_bins, len(genes))
    rank_sum, nnz, total = np.zeros(shape), np.zeros(shape), np.zeros(shape)
    tie = np.zeros(len(genes))
    for start in range(0, len(genes), chunk):
        block = matrix[:, genes[start:start + chunk]]
        block = sp.csc_matrix(block[cells] if cells is not None else block)
        for j in range(block.shape[1]):
            g = start + j
            lo, hi = block.indptr[j], block.indptr[j + 1]
            values, member = block.data[lo:hi].astype(np.float64), bins[block.indices[lo:hi]]
            uniq, inverse, counts = np.unique(np.append(values, 0.0), return_inverse=True, return_counts=True)
            zero = inverse[-1]
            counts[zero] += n - len(values) - 1  # the appended 0 stands for every implicit zero
            average_rank = np.cumsum(counts) - (counts - 1) / 2.0
            stored = np.bincount(member, minlength=n_bins)
            rank_sum[:, g] = np.bincount(member, weights=average_rank[inverse[:-1]], minlength=n_bins) + \
                (bin_counts - stored) * average_rank[zero]
            nnz[:, g] = np.bincount(member, weights=values != 0, minlength=n_bins)
            total[:, g] = np.bincount(member, weights=values, minlength=n_bins)
            tie[g] = (counts.astype(np.float64) ** 3 - counts).sum()
    n1 = n_group[:, None]
    n2 = n - n1
    with np.errstate(divide="ignore", invalid="ignore"):
        sd = np.sqrt(n1 * n2 / 12.0 * ((n + 1) - tie / (n * (n - 1.0))))
        z = (rank_sum[:k] - n1 * (n + 1) / 2.0) / sd
        mean_group = total[:k] / n1
        mean_rest = (total.sum(0) - total[:k]) / n2
        pct_group = nnz[:k] / n1
        pct_rest = (nnz.sum(0) - nnz[:k]) / n2
    z[~np.isfinite(z)] = 0.0  # constant among measured cells; scanpy also reports 0
    keep = np.broadcast_to((n1 >= max(min_measured, 1)) & (n2 >= max(min_measured, 1)), z.shape)
    group_index, gene_index = np.nonzero(keep)
    zk = z[keep]
    return pd.DataFrame({
        "group": categories[group_index].astype(str),
        "names": np.asarray(var_names)[genes][gene_index].astype(str),
        "scores": zk,
        "logfoldchanges": np.log2((np.expm1(mean_group[keep]) + 1e-9) / (np.expm1(mean_rest[keep]) + 1e-9)),
        "pct_nz_group": pct_group[keep],
        "pct_nz_reference": pct_rest[keep],
        "n_measured_group": n1[group_index, 0].astype(np.int64),
        "n_measured_rest": n2[group_index, 0].astype(np.int64),
        "U": rank_sum[:k][keep] - n1[group_index, 0] * (n1[group_index, 0] + 1) / 2.0,
        "log10_pval": (np.log(2.0) + norm.logsf(np.abs(zk))) / np.log(10.0),
    })


def _log10_bh(log10_p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg in log10 space, so q stays finite where p underflows."""
    m = log10_p.size
    order = np.argsort(log10_p)
    q = log10_p[order] + np.log10(m) - np.log10(np.arange(1, m + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 0.0)
    return out


def measured_everywhere_view(adata, measurement=None):
    """COSG-style scores have no per-gene cell mask: restrict them to genes every present sample measured.

    A partly measured gene would score its structural zeros as absence, so it
    is left out rather than scored. Returns ``adata`` (no measurement, or all
    genes measured) or a gene-subset view.
    """
    if measurement is None:
        return adata
    measurement.check(adata)
    full = measurement.measured_in_all(adata.obs[measurement.sample_key])
    return adata if full.all() else adata[:, full]


def run_markers_omicverse(adata, groupby: str = "cluster_identity", layer: str = "log1p_norm", method: str = "cosg", n_genes: int = 100,
                          measurement=None) -> pd.DataFrame:
    """Run OmicVerse marker wrapper and immediately export the result.

    OmicVerse stores marker results in adata.uns; this wrapper converts them
    to a DataFrame and removes the heavy intermediate key to follow the
    minimal-H5AD policy. With ``measurement``, only genes measured in every
    sample are scored (``measured_everywhere_view``).
    """
    from scsp_agent_sop.omicverse_facilities import require_omicverse

    ov = require_omicverse()
    view = measured_everywhere_view(adata, measurement)
    # OmicVerse writes uns; give it its own object rather than a view of the caller's
    adata = view.copy() if view.is_view else adata
    key = "_ov_markers_tmp"
    ov.single.find_markers(adata, groupby=groupby, method=method, n_genes=n_genes, key_added=key, layer=layer)
    # Prefer OmicVerse helper if present.
    try:
        df = ov.single.get_markers(adata, key=key)
    except Exception:
        import scanpy as sc
        df = sc.get.rank_genes_groups_df(adata, group=None, key=key)
    if key in adata.uns:
        del adata.uns[key]
    return df


def run_markers_fastde(adata, groupby: str = "cluster_identity", layer: str = "log1p_norm", method: str = "cosg", n_genes: int = 100,
                       measurement=None) -> pd.DataFrame:
    if method != "cosg":
        raise ValueError("FastDE annotation marker wrapper currently supports method='cosg'")
    if groupby not in adata.obs:
        raise KeyError(f"groupby {groupby!r} is not present in obs")
    from fastde.markers import run_cosg_markers

    adata = measured_everywhere_view(adata, measurement)

    matrix = adata.layers[layer] if layer in adata.layers else adata.X
    return run_cosg_markers(matrix, adata.obs[groupby], adata.var_names, top_n=n_genes).table


def build_annotation_evidence_template(adata, cluster_key: str = "cluster_identity") -> pd.DataFrame:
    rows = []
    for cl, idx in adata.obs.groupby(cluster_key).indices.items():
        obs = adata.obs.iloc[list(idx)]
        rows.append({
            "cluster": cl,
            "n_cells": len(obs),
            "n_samples": obs["sample_id"].nunique() if "sample_id" in obs else None,
            "proposed_label": "review_required",
            "marker_evidence": "pending",
            "nmf_evidence": "pending",
            "knn_evidence": "pending",
            "reference_evidence": "pending",
            "confidence": "low",
            "human_review_required": True,
        })
    return pd.DataFrame(rows)
