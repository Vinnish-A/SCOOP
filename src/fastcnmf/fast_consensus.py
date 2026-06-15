from __future__ import annotations

import multiprocessing as mp
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class ConsensusSummary:
    run_name: str
    output_dir: str
    k_values: tuple[int, ...]
    workers: int
    scheduler: str
    k_selection: bool


def _load_cnmf_symbols():
    from cnmf.cnmf import cNMF, load_df_from_npz, save_df_to_npz

    return cNMF, load_df_from_npz, save_df_to_npz


def _fast_prepare_manifest_path(output_dir: Path, run_name: str) -> Path:
    return output_dir / run_name / "cnmf_tmp" / f"{run_name}.fast_prepare_manifest.json"


def _read_names(path: Path) -> np.ndarray:
    return np.array(path.read_text(encoding="utf-8").splitlines(), dtype=object)


def _load_norm_counts_for_lite(obj, output_dir: Path, run_name: str):
    import json
    import scanpy as sc

    manifest_path = _fast_prepare_manifest_path(output_dir, run_name)
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        npy = manifest.get("normalized_counts_npy")
        obs_names = manifest.get("normalized_counts_obs_names")
        if npy and obs_names and Path(npy).is_file() and Path(obs_names).is_file():
            return np.load(npy, mmap_mode="r"), pd.Index(_read_names(Path(obs_names)))

    norm_counts = sc.read(obj.paths["normalized_counts"])
    return norm_counts.X, norm_counts.obs.index


def _load_matrix_cache(path: Path, fmt: str):
    if fmt == "sparse_npz":
        import scipy.sparse as sp

        return sp.load_npz(path)
    if fmt == "npy":
        return np.load(path, mmap_mode="r")
    raise ValueError(f"unsupported matrix cache format: {fmt}")


def _load_tpm_hvg_cache_for_lite(output_dir: Path, run_name: str):
    import json

    manifest_path = _fast_prepare_manifest_path(output_dir, run_name)
    if not manifest_path.is_file():
        return None
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    required = ["tpm_hvg_raw", "tpm_hvg_scaled", "tpm_hvg_var_names"]
    if not all(manifest.get(key) for key in required):
        return None
    raw_path = Path(manifest["tpm_hvg_raw"])
    scaled_path = Path(manifest["tpm_hvg_scaled"])
    var_path = Path(manifest["tpm_hvg_var_names"])
    if not raw_path.is_file() or not scaled_path.is_file() or not var_path.is_file():
        return None
    return {
        "raw": _load_matrix_cache(raw_path, manifest.get("tpm_hvg_raw_format", "sparse_npz")),
        "scaled": _load_matrix_cache(scaled_path, manifest.get("tpm_hvg_scaled_format", "sparse_npz")),
        "var_names": pd.Index(_read_names(var_path)),
    }


def _all_k_values(output_dir: Path, run_name: str) -> tuple[int, ...]:
    cNMF, load_df_from_npz, _ = _load_cnmf_symbols()
    obj = cNMF(output_dir=str(output_dir), name=run_name)
    params = load_df_from_npz(obj.paths["nmf_replicate_parameters"])
    return tuple(int(k) for k in sorted(set(params.n_components)))


def _run_consensus_k(job: dict) -> dict:
    cNMF, _, _ = _load_cnmf_symbols()
    obj = cNMF(output_dir=job["output_dir"], name=job["run_name"])
    k = int(job["k"])
    obj.combine(components=k)
    obj.consensus(
        k,
        density_threshold=job["density_threshold"],
        local_neighborhood_size=job["local_neighborhood_size"],
        show_clustering=job["show_clustering"],
        build_ref=job["build_reference"],
        close_clustergram_fig=True,
    )
    return {"k": k, "stage": "consensus"}


def _run_consensus_lite_k(job: dict) -> dict:
    import scanpy as sc
    import scipy.sparse as sp
    from sklearn.cluster import KMeans
    from sklearn.metrics.pairwise import euclidean_distances

    cNMF, load_df_from_npz, save_df_to_npz = _load_cnmf_symbols()
    obj = cNMF(output_dir=job["output_dir"], name=job["run_name"])
    k = int(job["k"])
    density_threshold = float(job["density_threshold"])
    density_threshold_repl = str(density_threshold).replace(".", "_")
    local_neighborhood_size = float(job["local_neighborhood_size"])

    obj.combine(components=k)
    merged_spectra = load_df_from_npz(obj.paths["merged_spectra"] % k)
    norm_counts_x, norm_counts_obs = _load_norm_counts_for_lite(
        obj,
        output_dir=Path(job["output_dir"]),
        run_name=job["run_name"],
    )
    n_neighbors = int(local_neighborhood_size * merged_spectra.shape[0] / k)
    if n_neighbors < 1:
        raise ValueError("lite consensus requires at least one local-density neighbor")

    l2_spectra = (merged_spectra.T / np.sqrt((merged_spectra**2).sum(axis=1))).T
    if density_threshold < 2:
        if Path(obj.paths["local_density_cache"] % k).is_file():
            local_density = load_df_from_npz(obj.paths["local_density_cache"] % k)
        else:
            topics_dist = euclidean_distances(l2_spectra.values)
            partitioning_order = np.argpartition(topics_dist, n_neighbors + 1)[:, : n_neighbors + 1]
            distance_to_nearest_neighbors = topics_dist[np.arange(topics_dist.shape[0])[:, None], partitioning_order]
            local_density = pd.DataFrame(
                distance_to_nearest_neighbors.sum(1) / n_neighbors,
                columns=["local_density"],
                index=l2_spectra.index,
            )
            save_df_to_npz(local_density, obj.paths["local_density_cache"] % k)
        density_filter = local_density.iloc[:, 0] < density_threshold
        l2_spectra = l2_spectra.loc[density_filter, :]
        if l2_spectra.shape[0] == 0:
            raise RuntimeError("Zero components remain after density filtering. Consider increasing density threshold")

    kmeans_model = KMeans(n_clusters=k, n_init=10, random_state=1)
    kmeans_model.fit(l2_spectra)
    kmeans_cluster_labels = pd.Series(kmeans_model.labels_ + 1, index=l2_spectra.index)

    median_spectra = l2_spectra.groupby(kmeans_cluster_labels).median()
    median_spectra = (median_spectra.T / median_spectra.sum(1)).T

    rf_usages = obj.refit_usage(norm_counts_x, median_spectra)
    rf_usages = pd.DataFrame(rf_usages, index=norm_counts_obs, columns=median_spectra.index)

    norm_usages = rf_usages.div(rf_usages.sum(axis=1), axis=0)
    reorder = norm_usages.sum(axis=0).sort_values(ascending=False)
    rf_usages = rf_usages.loc[:, reorder.index]
    norm_usages = norm_usages.loc[:, reorder.index]
    median_spectra = median_spectra.loc[reorder.index, :]
    rf_usages.columns = np.arange(1, rf_usages.shape[1] + 1)
    norm_usages.columns = rf_usages.columns
    median_spectra.index = rf_usages.columns

    tpm_stats = load_df_from_npz(obj.paths["tpm_stats"])
    hvgs = open(obj.paths["nmf_genes_list"], encoding="utf-8").read().split("\n")

    tpm_hvg_cache = _load_tpm_hvg_cache_for_lite(Path(job["output_dir"]), job["run_name"])
    if tpm_hvg_cache is not None:
        raw_tpm_hvg = tpm_hvg_cache["raw"]
        norm_tpm_x = tpm_hvg_cache["scaled"]
        tpm_hvg_var = tpm_hvg_cache["var_names"]
    else:
        tpm = sc.read(obj.paths["tpm"])
        tpm_hvg = tpm[:, hvgs]
        raw_tpm_hvg = tpm_hvg.X
        tpm_hvg_var = pd.Index(tpm_hvg.var.index)
        norm_tpm = tpm_hvg.copy()
        if sp.issparse(norm_tpm.X):
            sc.pp.scale(norm_tpm, zero_center=False)
        else:
            norm_tpm.X /= norm_tpm.X.std(axis=0, ddof=1)
        norm_tpm_x = norm_tpm.X

    spectra_tpm_rf = obj.refit_spectra(raw_tpm_hvg, norm_usages.astype(raw_tpm_hvg.dtype))
    spectra_tpm_rf = pd.DataFrame(spectra_tpm_rf, index=rf_usages.columns, columns=tpm_hvg_var)
    spectra_tpm_rf = spectra_tpm_rf.div(tpm_stats.loc[hvgs, "__std"], axis=1)
    rf_usages = obj.refit_usage(norm_tpm_x, spectra_tpm_rf.astype(norm_tpm_x.dtype))
    rf_usages = pd.DataFrame(rf_usages, index=norm_counts_obs, columns=spectra_tpm_rf.index)

    save_df_to_npz(median_spectra, obj.paths["consensus_spectra"] % (k, density_threshold_repl))
    save_df_to_npz(rf_usages, obj.paths["consensus_usages"] % (k, density_threshold_repl))
    median_spectra.to_csv(obj.paths["consensus_spectra__txt"] % (k, density_threshold_repl), sep="\t")
    rf_usages.to_csv(obj.paths["consensus_usages__txt"] % (k, density_threshold_repl), sep="\t")
    return {"k": k, "stage": "consensus-lite"}


def _run_k_selection_stats(job: dict) -> tuple[int, pd.Series]:
    cNMF, _, _ = _load_cnmf_symbols()
    obj = cNMF(output_dir=job["output_dir"], name=job["run_name"])
    stats = obj.consensus(
        int(job["k"]),
        density_threshold=job["density_threshold"],
        local_neighborhood_size=job["local_neighborhood_size"],
        show_clustering=False,
        build_ref=False,
        skip_density_and_return_after_stats=True,
        close_clustergram_fig=True,
    )
    return int(job["k"]), stats.stats


def _write_k_selection_plot(output_dir: Path, run_name: str, stats: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt

    cNMF, _, save_df_to_npz = _load_cnmf_symbols()
    obj = cNMF(output_dir=str(output_dir), name=run_name)
    stats = stats.sort_values("k").reset_index(drop=True)
    save_df_to_npz(stats, obj.paths["k_selection_stats"])

    fig = plt.figure(figsize=(6, 4))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    ax1.plot(stats.k, stats.silhouette, "o-", color="b")
    ax1.set_ylabel("Stability", color="b", fontsize=15)
    for tick in ax1.get_yticklabels():
        tick.set_color("b")
    ax2.plot(stats.k, stats.prediction_error, "o-", color="r")
    ax2.set_ylabel("Error", color="r", fontsize=15)
    for tick in ax2.get_yticklabels():
        tick.set_color("r")
    ax1.set_xlabel("Number of Components", fontsize=15)
    ax1.grid("on")
    plt.tight_layout()
    fig.savefig(obj.paths["k_selection_plot"], dpi=250)
    plt.close(fig)


def _run_jobs(jobs: list[dict], workers: int, func) -> list:
    if workers == 1:
        return [func(job) for job in jobs]
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=workers) as pool:
        return list(pool.imap_unordered(func, jobs, chunksize=1))


def run_fast_consensus(
    *,
    output_dir: Path,
    run_name: str,
    workers: int,
    k_values: tuple[int, ...] | None = None,
    density_threshold: float = 0.5,
    local_neighborhood_size: float = 0.30,
    show_clustering: bool = False,
    build_reference: bool = True,
    k_selection: bool = True,
    lite: bool = False,
) -> dict:
    """Run cNMF-compatible consensus with k-level parallel scheduling."""

    if workers < 1:
        raise ValueError("workers must be positive")
    if k_values is None or len(k_values) == 0:
        k_values = _all_k_values(output_dir, run_name)
    worker_count = min(workers, len(k_values))
    base_job = {
        "output_dir": str(output_dir),
        "run_name": run_name,
        "density_threshold": density_threshold,
        "local_neighborhood_size": local_neighborhood_size,
        "show_clustering": show_clustering,
        "build_reference": build_reference,
    }
    jobs = [{**base_job, "k": int(k)} for k in k_values]
    consensus_results = _run_jobs(jobs, worker_count, _run_consensus_lite_k if lite else _run_consensus_k)

    k_selection_results = []
    if k_selection and not lite:
        stats_pairs = _run_jobs(jobs, worker_count, _run_k_selection_stats)
        stats = pd.DataFrame([series for _, series in sorted(stats_pairs, key=lambda item: item[0])])
        _write_k_selection_plot(output_dir=output_dir, run_name=run_name, stats=stats)
        k_selection_results = [{"k": int(k), "stage": "k_selection"} for k, _ in stats_pairs]

    summary = ConsensusSummary(
        run_name=run_name,
        output_dir=str(output_dir),
        k_values=tuple(int(k) for k in k_values),
        workers=worker_count,
        scheduler="parallel-by-k-lite" if lite else "parallel-by-k",
        k_selection=k_selection and not lite,
    )
    return {
        **summary.__dict__,
        "lite": lite,
        "consensus_results": sorted(consensus_results, key=lambda item: item["k"]),
        "k_selection_results": sorted(k_selection_results, key=lambda item: item["k"]),
    }
