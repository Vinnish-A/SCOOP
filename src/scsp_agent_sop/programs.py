from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import sparse

from .matrix import get_matrix


FASTCNMF_DEFAULT_MAX_ITER = 50
FASTCNMF_DEFAULT_SEEDS = tuple(range(20))
FASTCNMF_COMPATIBLE_METHODS = {"fastcnmf", "fast_consensus_nmf", "sklearn_nmf"}


@dataclass
class NMFRunResult:
    k: int
    seed: int
    error: float
    W: np.ndarray
    H: np.ndarray


def _as_nonnegative_matrix(X):
    if sparse.issparse(X):
        if X.min() < 0:
            raise ValueError("NMF input contains negative values")
        return X
    X = np.asarray(X, dtype=float)
    if np.nanmin(X) < 0:
        raise ValueError("NMF input contains negative values")
    return X


def _normalize_components(H: np.ndarray) -> np.ndarray:
    denom = np.linalg.norm(H, axis=1, keepdims=True)
    return H / np.maximum(denom, 1e-12)


def run_nmf_once(X, k: int, seed: int, max_iter: int = FASTCNMF_DEFAULT_MAX_ITER) -> NMFRunResult:
    from sklearn.decomposition import NMF

    model = NMF(
        n_components=k,
        init="nndsvda",
        solver="cd",
        beta_loss="frobenius",
        max_iter=max_iter,
        tol=1e-4,
        random_state=seed,
    )
    W = model.fit_transform(X)
    H = model.components_
    return NMFRunResult(k=k, seed=seed, error=float(model.reconstruction_err_), W=W, H=H)


def consensus_stability(results: list[NMFRunResult], k: int) -> tuple[float, np.ndarray]:
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.metrics.pairwise import cosine_similarity

    H_all = np.vstack([_normalize_components(r.H) for r in results if r.k == k])
    if H_all.shape[0] <= k:
        return 1.0, H_all
    sim = cosine_similarity(H_all)
    dist = 1 - np.clip(sim, 0, 1)
    clustering = AgglomerativeClustering(n_clusters=k, metric="precomputed", linkage="average")
    labels = clustering.fit_predict(dist)
    stabilities = []
    consensus = []
    for label in range(k):
        idx = np.where(labels == label)[0]
        if len(idx) == 0:
            continue
        sub = sim[np.ix_(idx, idx)]
        stabilities.append(float(np.mean(sub[np.triu_indices_from(sub, k=1)])) if len(idx) > 1 else 1.0)
        consensus.append(H_all[idx].mean(axis=0))
    return float(np.median(stabilities)) if stabilities else 0.0, np.vstack(consensus)


def choose_k(summary: pd.DataFrame, stability_threshold: float = 0.70) -> int:
    stable = summary[summary["median_stability"] >= stability_threshold].copy()
    if stable.empty:
        return int(summary.sort_values(["median_stability", "median_error"], ascending=[False, True]).iloc[0]["k"])
    # Prefer the smallest stable K after error has begun flattening.
    stable = stable.sort_values("k")
    return int(stable.iloc[0]["k"])


def run_fast_consensus_nmf(
    adata,
    *,
    layer: str = "log1p_norm",
    hvg_key: str = "highly_variable_biology",
    k_grid: Iterable[int] = (5, 8, 10, 12, 15),
    seeds: Iterable[int] = FASTCNMF_DEFAULT_SEEDS,
    max_iter: int = FASTCNMF_DEFAULT_MAX_ITER,
    stability_threshold: float = 0.70,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if hvg_key in adata.var:
        gene_mask = adata.var[hvg_key].to_numpy(bool)
    else:
        gene_mask = np.ones(adata.n_vars, dtype=bool)
    # Exclude genes that should not drive programme discovery by default.
    for col in ["mt_gene", "hb_gene"]:
        if col in adata.var:
            gene_mask &= ~adata.var[col].to_numpy(bool)
    genes = np.asarray(adata.var_names)[gene_mask]
    X = _as_nonnegative_matrix(get_matrix(adata, layer)[:, gene_mask])
    results: list[NMFRunResult] = []
    for k in k_grid:
        for seed in seeds:
            results.append(run_nmf_once(X, int(k), int(seed), max_iter=max_iter))
    summary_rows = []
    for k in k_grid:
        kr = [r for r in results if r.k == int(k)]
        stability, _ = consensus_stability(results, int(k))
        summary_rows.append({
            "k": int(k),
            "median_error": float(np.median([r.error for r in kr])),
            "min_error": float(np.min([r.error for r in kr])),
            "median_stability": stability,
            "n_seeds": len(kr),
        })
    summary = pd.DataFrame(summary_rows)
    selected_k = choose_k(summary, stability_threshold=stability_threshold)
    best = min([r for r in results if r.k == selected_k], key=lambda r: r.error)
    adata.obsm["X_nmf_usage"] = best.W
    adata.obs["dominant_nmf_program"] = pd.Categorical([f"P{int(i)+1}" for i in np.argmax(best.W, axis=1)])
    usage_sum = np.maximum(best.W.sum(axis=1), 1e-12)
    probs = best.W / usage_sum[:, None]
    adata.obs["nmf_program_entropy"] = -(probs * np.log(np.maximum(probs, 1e-12))).sum(axis=1)
    weights = pd.DataFrame(best.H, columns=genes)
    weights.insert(0, "program", [f"P{i+1}" for i in range(best.H.shape[0])])
    usage = pd.DataFrame(best.W, index=adata.obs_names, columns=[f"P{i+1}" for i in range(best.W.shape[1])])
    usage.insert(0, "obs_name", usage.index)
    summary["selected"] = summary["k"] == selected_k
    return summary, weights, usage


def run_omicverse_cnmf_validation(adata, components, output_dir: str, n_iter: int = 100, use_gpu: bool = True):
    """Optional validation path using OmicVerse cNMF.

    This is intentionally not the default. It is called only when the fast
    consensus NMF is unstable or the programme supports an important claim.
    """
    import omicverse as ov
    return ov.single.cNMF(
        adata,
        components=components,
        n_iter=n_iter,
        output_dir=output_dir,
        use_gpu=use_gpu,
        name="omicverse_cnmf_validation",
    )
