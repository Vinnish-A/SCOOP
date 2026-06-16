from __future__ import annotations

from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd


def normalize_log1p(adata, counts_layer: str = "counts", target_sum: float = 1e4, layer_out: str = "log1p_norm") -> None:
    import scanpy as sc
    adata.X = adata.layers[counts_layer].copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.layers[layer_out] = adata.X.copy()


def select_hvg(
    adata,
    counts_layer: str = "counts",
    batch_key: str | None = "sample_id",
    flavor: str = "seurat_v3",
    n_top_genes: int = 3000,
    output_key: str = "highly_variable_biology",
) -> pd.DataFrame:
    import scanpy as sc
    try:
        sc.pp.highly_variable_genes(
            adata,
            layer=counts_layer,
            batch_key=batch_key if batch_key in adata.obs else None,
            flavor=flavor,
            n_top_genes=n_top_genes,
            inplace=True,
        )
    except Exception:
        # seurat_v3 requires scikit-misc in some environments. Fallback is still deterministic.
        sc.pp.highly_variable_genes(
            adata,
            layer=None,
            batch_key=batch_key if batch_key in adata.obs else None,
            flavor="cell_ranger",
            n_top_genes=n_top_genes,
            inplace=True,
        )
        flavor = "cell_ranger_fallback"
    adata.var[output_key] = adata.var["highly_variable"].astype(bool)
    rank_cols = [c for c in ["highly_variable", "highly_variable_rank", "means", "variances", "variances_norm"] if c in adata.var]
    out = adata.var[rank_cols].copy()
    out.insert(0, "gene", adata.var_names)
    out["hvg_key"] = output_key
    out["method"] = flavor
    return out


def run_pca(adata, hvg_key: str, obsm_key: str, n_comps: int = 50, layer: str = "log1p_norm") -> None:
    import scanpy as sc
    mask = adata.var[hvg_key].to_numpy(bool)
    ad = adata[:, mask].copy()
    ad.X = ad.layers[layer].copy() if layer in ad.layers else ad.X.copy()
    sc.pp.pca(ad, n_comps=n_comps, svd_solver="arpack")
    adata.obsm[obsm_key] = ad.obsm["X_pca"].copy()
    # Keep loadings only as small matrix if needed by users; default is external diagnostics.


def score_programs(adata, organism: str = "human", layer: str = "log1p_norm") -> None:
    import scanpy as sc
    old = adata.X
    if layer in adata.layers:
        adata.X = adata.layers[layer].copy()
    stress = ["FOS", "JUN", "JUNB", "JUND", "ATF3", "EGR1", "DUSP1", "HSPA1A", "HSPA1B", "DNAJB1"]
    cycle_s = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"]
    cycle_g2m = ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"]
    if organism.lower().startswith("mouse"):
        def mm(genes): return [g.capitalize() for g in genes]
        stress, cycle_s, cycle_g2m = map(mm, [stress, cycle_s, cycle_g2m])
    genes = set(map(str, adata.var_names))
    sc.tl.score_genes(adata, [g for g in stress if g in genes], score_name="stress_score")
    ribo_genes = list(adata.var_names[adata.var.get("ribo_gene", False).to_numpy(bool)]) if "ribo_gene" in adata.var else []
    if ribo_genes:
        sc.tl.score_genes(adata, ribo_genes, score_name="ribo_score")
    else:
        adata.obs["ribo_score"] = 0.0
    sc.tl.score_genes_cell_cycle(
        adata,
        s_genes=[g for g in cycle_s if g in genes],
        g2m_genes=[g for g in cycle_g2m if g in genes],
    )
    adata.obs.rename(columns={"S_score": "cell_cycle_s_score", "G2M_score": "cell_cycle_g2m_score"}, inplace=True)
    adata.obs["proliferation_score"] = adata.obs[["cell_cycle_s_score", "cell_cycle_g2m_score"]].max(axis=1)
    adata.X = old


def robust_z_by_sample(adata, cols: Sequence[str], sample_key: str = "sample_id") -> None:
    for col in cols:
        out = np.zeros(adata.n_obs, dtype=float)
        for _, idx in adata.obs.groupby(sample_key).indices.items():
            x = adata.obs.iloc[list(idx)][col].to_numpy(float)
            med = np.nanmedian(x)
            mad = max(1.4826 * np.nanmedian(np.abs(x - med)), 1e-9)
            out[list(idx)] = (x - med) / mad
        adata.obs[f"{col}_z_sample"] = out


def build_identity_hvg_from_program_decision(adata, biology_key: str = "highly_variable_biology", output_key: str = "highly_variable_identity") -> pd.DataFrame:
    mask = adata.var[biology_key].to_numpy(bool).copy()
    excluded = np.zeros(adata.n_vars, dtype=bool)
    # Minimal rule: always exclude mt genes from identity HVG. Other programme genes are
    # excluded only if earlier diagnostics set exclude_from_identity_hvg.
    for col in ["mt_gene", "exclude_from_identity_hvg"]:
        if col in adata.var:
            excluded |= adata.var[col].to_numpy(bool)
    mask &= ~excluded
    adata.var[output_key] = mask
    df = pd.DataFrame({"gene": adata.var_names, output_key: mask, "excluded_from_identity": excluded})
    return df


def run_harmony2(
    adata,
    basis: str = "X_pca_identity_prebatch",
    batch_keys: str | list[str] = "sample_id",
    output: str = "X_pca_harmony_identity",
    max_iter_harmony: int = 20,
    random_state: int = 0,
    ncores: int = 0,
    sigma: float = 0.1,
) -> None:
    try:
        import harmonypy
    except Exception as exc:
        raise ImportError("Harmony 2.0 is not available. Install harmonypy>=2.0,<3.") from exc
    keys = [batch_keys] if isinstance(batch_keys, str) else list(batch_keys)
    keys = [k for k in keys if k in adata.obs and adata.obs[k].nunique() > 1]
    if not keys:
        adata.obsm[output] = adata.obsm[basis].copy()
        return
    nclust = int(min(round(adata.n_obs / 30.0), 100))
    nclust = max(nclust, 1)
    result = harmonypy.run_harmony(
        np.asarray(adata.obsm[basis], dtype=np.float64),
        adata.obs,
        keys,
        sigma=np.repeat(float(sigma), nclust).astype(np.float64),
        nclust=nclust,
        max_iter_harmony=int(max_iter_harmony),
        verbose=False,
        random_state=int(random_state),
        ncores=int(ncores),
    )
    corrected = np.asarray(result.Z_corr, dtype=np.float64)
    if corrected.shape != adata.obsm[basis].shape and corrected.T.shape == adata.obsm[basis].shape:
        corrected = corrected.T
    if corrected.shape != adata.obsm[basis].shape:
        raise ValueError(f"Harmony 2.0 returned shape {corrected.shape}, expected {adata.obsm[basis].shape}.")
    adata.obsm[output] = corrected

def neighbors_umap(adata, use_rep: str, prefix: str, n_neighbors: int = 15, n_pcs: int | None = None, min_dist: float = 0.3, random_state: int = 0) -> None:
    import scanpy as sc
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors, n_pcs=n_pcs, key_added=f"neighbors_{prefix}")
    # Copy graph to stable names.
    adata.obsp[f"connectivities_{prefix}"] = adata.obsp[f"neighbors_{prefix}_connectivities"].copy()
    adata.obsp[f"distances_{prefix}"] = adata.obsp[f"neighbors_{prefix}_distances"].copy()
    sc.tl.umap(adata, neighbors_key=f"neighbors_{prefix}", min_dist=min_dist, random_state=random_state)
    adata.obsm[f"X_umap_{prefix}"] = adata.obsm["X_umap"].copy()


def _unique_sorted(values: Iterable[float]) -> list[float]:
    return sorted({round(float(value), 6) for value in values})


def _resolution_neighborhood(center: float, *, lower: float, upper: float, window: float, step: float) -> list[float]:
    start = max(lower, center - window)
    stop = min(upper, center + window)
    n_steps = int(np.floor((stop - start) / step + 0.5))
    values = [start + i * step for i in range(n_steps + 1)]
    values.append(center)
    return _unique_sorted(value for value in values if lower <= value <= upper)


def _choose_resolution(stability: pd.DataFrame, *, min_ari: float, min_clusters: int, max_clusters: int) -> float:
    candidates = stability[
        (stability["median_ari"] >= min_ari)
        & (stability["n_clusters_seed0"] >= min_clusters)
        & (stability["n_clusters_seed0"] <= max_clusters)
    ]
    if len(candidates):
        return float(candidates.sort_values(["resolution"]).iloc[0]["resolution"])
    ranked = stability.sort_values(["median_ari", "n_clusters_seed0", "resolution"], ascending=[False, False, True])
    return float(ranked.iloc[0]["resolution"])


def leiden_sweep(
    adata,
    graph_prefix: str = "identity",
    resolutions: Iterable[float] = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5),
    seeds: Iterable[int] = (0, 1, 2, 3, 4),
    search_config: Mapping[str, Any] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    import scanpy as sc
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    cfg = dict(search_config or {})
    strategy = str(cfg.get("strategy", "coarse_to_fine"))
    all_resolutions = _unique_sorted(resolutions)
    all_seeds = [int(seed) for seed in seeds]
    min_ari = float(cfg.get("min_seed_ari", 0.8))
    min_clusters = int(cfg.get("min_clusters", 1))
    max_clusters = int(cfg.get("max_clusters", 10**9))

    sweep_rows = []
    labels_by_res: dict[float, list[np.ndarray]] = {}
    seed_by_res: dict[float, list[int]] = {}

    def run_jobs(job_resolutions: Iterable[float], job_seeds: Iterable[int], phase: str) -> None:
        for res in _unique_sorted(job_resolutions):
            labels_by_res.setdefault(float(res), [])
            seed_by_res.setdefault(float(res), [])
            existing = set(seed_by_res[float(res)])
            for seed in [int(value) for value in job_seeds]:
                if seed in existing:
                    continue
                key = f"_tmp_leiden_r{res}_s{seed}"
                sc.tl.leiden(
                    adata,
                    resolution=float(res),
                    random_state=int(seed),
                    key_added=key,
                    adjacency=adata.obsp[f"connectivities_{graph_prefix}"],
                )
                labels = adata.obs[key].astype(str).to_numpy()
                labels_by_res[float(res)].append(labels)
                seed_by_res[float(res)].append(int(seed))
                sweep_rows.append({
                    "resolution": float(res),
                    "seed": int(seed),
                    "n_clusters": int(pd.Series(labels).nunique()),
                    "phase": phase,
                })

    def summarize() -> pd.DataFrame:
        stability_rows = []
        for res, labels_list in labels_by_res.items():
            ari = []
            nmi = []
            for i in range(len(labels_list)):
                for j in range(i + 1, len(labels_list)):
                    ari.append(adjusted_rand_score(labels_list[i], labels_list[j]))
                    nmi.append(normalized_mutual_info_score(labels_list[i], labels_list[j]))
            stability_rows.append({
                "resolution": res,
                "median_ari": float(np.median(ari)) if ari else 1.0,
                "median_nmi": float(np.median(nmi)) if nmi else 1.0,
                "n_clusters_seed0": int(pd.Series(labels_by_res[res][0]).nunique()),
                "n_seeds": int(len(labels_list)),
                "search_strategy": strategy,
            })
        return pd.DataFrame(stability_rows).sort_values(["median_ari", "resolution"], ascending=[False, True])

    if strategy == "coarse_to_fine" and len(all_resolutions) > 2:
        coarse_resolutions = cfg.get("coarse_resolutions", [all_resolutions[0], all_resolutions[len(all_resolutions) // 2], all_resolutions[-1]])
        coarse_resolutions = [float(res) for res in coarse_resolutions if min(all_resolutions) <= float(res) <= max(all_resolutions)]
        coarse_seeds = [int(seed) for seed in cfg.get("coarse_seeds", all_seeds[:2] or [0])]
        run_jobs(coarse_resolutions, coarse_seeds, "coarse")
        coarse_stability = summarize()
        center = _choose_resolution(coarse_stability, min_ari=min_ari, min_clusters=min_clusters, max_clusters=max_clusters)
        refine_resolutions = _resolution_neighborhood(
            center,
            lower=min(all_resolutions),
            upper=max(all_resolutions),
            window=float(cfg.get("refine_window", 0.25)),
            step=float(cfg.get("refine_step", 0.125)),
        )
        run_jobs(refine_resolutions, all_seeds, "refine")
    else:
        run_jobs(all_resolutions, all_seeds, "full")

    stability = summarize()
    full_seed_count = max(1, len(all_seeds))
    final_pool = stability[stability["n_seeds"] >= min(full_seed_count, max(1, len(all_seeds)))]
    if len(final_pool) == 0:
        final_pool = stability
    chosen = _choose_resolution(final_pool, min_ari=min_ari, min_clusters=min_clusters, max_clusters=max_clusters)
    chosen_key = f"_tmp_leiden_r{chosen}_s0"
    if chosen_key not in adata.obs:
        fallback_seed = seed_by_res[chosen][0]
        chosen_key = f"_tmp_leiden_r{chosen}_s{fallback_seed}"
    adata.obs["cluster_identity"] = adata.obs[chosen_key].astype("category")
    # Clean tmp columns from obs; keep full sweep externally.
    for c in list(adata.obs.columns):
        if c.startswith("_tmp_leiden_"):
            del adata.obs[c]
    sweep = pd.DataFrame(sweep_rows)
    stability["chosen"] = stability["resolution"] == chosen
    return sweep, stability
