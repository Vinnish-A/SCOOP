from __future__ import annotations

from pathlib import Path
import pandas as pd


def run_markers_scanpy(adata, groupby: str = "cluster_identity", layer: str = "log1p_norm", method: str = "wilcoxon", n_genes: int = 100) -> pd.DataFrame:
    import scanpy as sc
    old = adata.X
    if layer in adata.layers:
        adata.X = adata.layers[layer].copy()
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, n_genes=n_genes, key_added="_markers_tmp")
    df = sc.get.rank_genes_groups_df(adata, group=None, key="_markers_tmp")
    if "_markers_tmp" in adata.uns:
        del adata.uns["_markers_tmp"]
    adata.X = old
    return df


def run_markers_omicverse(adata, groupby: str = "cluster_identity", layer: str = "log1p_norm", method: str = "cosg", n_genes: int = 100) -> pd.DataFrame:
    """Run OmicVerse marker wrapper and immediately export the result.

    OmicVerse stores marker results in adata.uns; this wrapper converts them
    to a DataFrame and removes the heavy intermediate key to follow the
    minimal-H5AD policy.
    """
    from scsp_agent_sop.omicverse_facilities import require_omicverse

    ov = require_omicverse()
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


def run_markers_fastde(adata, groupby: str = "cluster_identity", layer: str = "log1p_norm", method: str = "cosg", n_genes: int = 100) -> pd.DataFrame:
    if method != "cosg":
        raise ValueError("FastDE annotation marker wrapper currently supports method='cosg'")
    if groupby not in adata.obs:
        raise KeyError(f"groupby {groupby!r} is not present in obs")
    from fastde.markers import run_cosg_markers

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
