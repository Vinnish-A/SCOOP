from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def _save_df_to_npz(obj: pd.DataFrame, path: Path, *, compressed: bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = np.savez_compressed if compressed else np.savez
    writer(path, data=obj.values, index=obj.index.values, columns=obj.columns.values)


def _load_df_from_npz(path: Path) -> pd.DataFrame:
    with np.load(path, allow_pickle=True) as f:
        return pd.DataFrame(**f)


def _link_or_copy(src: Path, dest: Path) -> str:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    try:
        os.link(src, dest)
        return "hardlink"
    except OSError:
        try:
            os.symlink(src.resolve(), dest)
            return "symlink"
        except OSError:
            shutil.copy2(src, dest)
            return "copy"


def _scale_zero_center_false(x: np.ndarray, *, output_dtype: str = "float64") -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    std = x.std(axis=0, ddof=1)
    std[std == 0] = 1.0
    x /= std
    if np.isnan(x).any():
        raise ValueError("NaNs produced while scaling normalized counts")
    if output_dtype == "float32":
        x = x.astype(np.float32)
    elif output_dtype != "float64":
        raise ValueError("output_dtype must be float32 or float64")
    return x


def _sparse_mean_var(matrix) -> tuple[np.ndarray, np.ndarray]:
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler(with_mean=False)
    scaler.fit(matrix)
    return np.asarray(scaler.mean_), np.asarray(scaler.var_)


def _dense_mean_var(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(matrix, dtype=np.float64)
    return x.mean(axis=0), x.var(axis=0, ddof=0)


def _nmf_iter_params(k_values: tuple[int, ...], n_iter: int, seed: int, output_template: str) -> pd.DataFrame:
    import itertools

    k_list = sorted(set(k_values))
    np.random.seed(seed=seed)
    seeds = np.random.randint(low=1, high=(2**31) - 1, size=len(k_values) * n_iter)
    rows = []
    for idx, (k, iteration) in enumerate(itertools.product(k_list, range(n_iter))):
        completed = Path(output_template % (k, iteration)).exists()
        rows.append([k, iteration, int(seeds[idx]), completed])
    return pd.DataFrame(rows, columns=["n_components", "iter", "nmf_seed", "completed"])


def run_fast_prepare(
    *,
    corrected_h5ad: Path,
    tpm_h5ad: Path | None,
    hvgs_txt: Path,
    output_dir: Path,
    run_name: str,
    k_values: tuple[int, ...],
    n_iter: int,
    seed: int = 20260614,
    beta_loss: str = "frobenius",
    max_nmf_iter: int = 50,
    init: str = "random",
    alpha_usage: float = 0.0,
    alpha_spectra: float = 0.0,
    norm_dtype: str = "float64",
    norm_store: str = "h5ad",
    precomputed_norm_npy: Path | None = None,
    precomputed_norm_obs_names: Path | None = None,
    precomputed_norm_var_names: Path | None = None,
    precomputed_tpm_stats: Path | None = None,
    precomputed_tpm_hvg_raw: Path | None = None,
    precomputed_tpm_hvg_scaled: Path | None = None,
    precomputed_tpm_hvg_obs_names: Path | None = None,
    precomputed_tpm_hvg_var_names: Path | None = None,
) -> dict:
    """Write cNMF-compatible prepare artifacts without invoking cNMF prepare."""

    import scanpy as sc
    from scipy import sparse

    if norm_store not in {"h5ad", "npy", "both"}:
        raise ValueError("norm_store must be h5ad, npy, or both")

    run_dir = output_dir / run_name
    tmp_dir = run_dir / "cnmf_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    hvgs = hvgs_txt.read_text(encoding="utf-8").rstrip().split("\n")
    use_precomputed_norm = precomputed_norm_npy is not None
    if use_precomputed_norm:
        if norm_store == "h5ad":
            raise ValueError("precomputed norm arrays require --norm-store npy or both")
        if precomputed_norm_obs_names is None or precomputed_norm_var_names is None:
            raise ValueError("precomputed norm arrays require obs and var name files")
        norm_x = np.load(precomputed_norm_npy, mmap_mode="r")
        norm_obs_names = precomputed_norm_obs_names.read_text(encoding="utf-8").splitlines()
        norm_var_names = precomputed_norm_var_names.read_text(encoding="utf-8").splitlines()
    else:
        corrected = sc.read_h5ad(corrected_h5ad)
        if corrected.n_vars != len(hvgs) or list(corrected.var_names) != hvgs:
            corrected = corrected[:, hvgs].copy()
        if sparse.issparse(corrected.X):
            norm_x = corrected.X.toarray()
        else:
            norm_x = np.asarray(corrected.X, dtype=np.float64).copy()
        norm_x = _scale_zero_center_false(norm_x, output_dtype=norm_dtype)
        norm_obs_names = list(map(str, corrected.obs_names))
        norm_var_names = list(map(str, corrected.var_names))

    tpm_n_vars = None
    if precomputed_tpm_stats is None:
        if tpm_h5ad is None:
            raise ValueError("tpm_h5ad is required when precomputed_tpm_stats is not provided")
        tpm = sc.read_h5ad(tpm_h5ad)
        tpm_n_vars = int(tpm.n_vars)
        if sparse.issparse(tpm.X):
            mean, var = _sparse_mean_var(tpm.X)
        else:
            mean, var = _dense_mean_var(tpm.X)
        tpm_stats = pd.DataFrame(
            [mean, np.sqrt(var)],
            index=["__mean", "__std"],
            columns=tpm.var.index,
        ).T

    normalized_counts_path = tmp_dir / f"{run_name}.norm_counts.h5ad"
    normalized_counts_npy_path = tmp_dir / f"{run_name}.norm_counts.npy"
    normalized_counts_obs_path = tmp_dir / f"{run_name}.norm_counts.obs.txt"
    normalized_counts_var_path = tmp_dir / f"{run_name}.norm_counts.var.txt"
    tpm_path = tmp_dir / f"{run_name}.tpm.h5ad"
    tpm_hvg_raw_path = tmp_dir / f"{run_name}.tpm_hvg.raw{precomputed_tpm_hvg_raw.suffix if precomputed_tpm_hvg_raw else '.npz'}"
    tpm_hvg_scaled_path = tmp_dir / f"{run_name}.tpm_hvg.scaled{precomputed_tpm_hvg_scaled.suffix if precomputed_tpm_hvg_scaled else '.npz'}"
    tpm_hvg_obs_path = tmp_dir / f"{run_name}.tpm_hvg.obs.txt"
    tpm_hvg_var_path = tmp_dir / f"{run_name}.tpm_hvg.var.txt"
    tpm_stats_path = tmp_dir / f"{run_name}.tpm_stats.df.npz"
    nmf_params_path = tmp_dir / f"{run_name}.nmf_params.df.npz"
    nmf_run_params_path = tmp_dir / f"{run_name}.nmf_idvrun_params.yaml"
    hvgs_out = run_dir / f"{run_name}.overdispersed_genes.txt"
    iter_template = str(tmp_dir / f"{run_name}.spectra.k_%d.iter_%d.df.npz")

    h5ad_written = False
    npy_written = False
    if norm_store in {"h5ad", "both"}:
        import anndata as ad

        norm_counts = ad.AnnData(
            X=np.asarray(norm_x),
            obs=pd.DataFrame(index=norm_obs_names),
            var=pd.DataFrame(index=norm_var_names),
        )
        norm_counts.write_h5ad(normalized_counts_path)
        h5ad_written = True
    if norm_store in {"npy", "both"}:
        if use_precomputed_norm:
            _link_or_copy(precomputed_norm_npy, normalized_counts_npy_path)
            _link_or_copy(precomputed_norm_obs_names, normalized_counts_obs_path)
            _link_or_copy(precomputed_norm_var_names, normalized_counts_var_path)
        else:
            np.save(normalized_counts_npy_path, norm_x)
            normalized_counts_obs_path.write_text("\n".join(norm_obs_names), encoding="utf-8")
            normalized_counts_var_path.write_text("\n".join(norm_var_names), encoding="utf-8")
        npy_written = True
    if tpm_h5ad is not None:
        tpm_link_mode = _link_or_copy(tpm_h5ad, tpm_path)
        tpm_out = str(tpm_path)
    else:
        tpm_link_mode = "none"
        tpm_out = None
    tpm_hvg_cache = {}
    if precomputed_tpm_hvg_raw is not None or precomputed_tpm_hvg_scaled is not None:
        if (
            precomputed_tpm_hvg_raw is None
            or precomputed_tpm_hvg_scaled is None
            or precomputed_tpm_hvg_obs_names is None
            or precomputed_tpm_hvg_var_names is None
        ):
            raise ValueError("precomputed TP10K HVG cache requires raw, scaled, obs, and var files")
        _link_or_copy(precomputed_tpm_hvg_raw, tpm_hvg_raw_path)
        _link_or_copy(precomputed_tpm_hvg_scaled, tpm_hvg_scaled_path)
        _link_or_copy(precomputed_tpm_hvg_obs_names, tpm_hvg_obs_path)
        _link_or_copy(precomputed_tpm_hvg_var_names, tpm_hvg_var_path)
        tpm_hvg_cache = {
            "tpm_hvg_raw": str(tpm_hvg_raw_path),
            "tpm_hvg_raw_format": "sparse_npz" if tpm_hvg_raw_path.suffix == ".npz" else "npy",
            "tpm_hvg_scaled": str(tpm_hvg_scaled_path),
            "tpm_hvg_scaled_format": "sparse_npz" if tpm_hvg_scaled_path.suffix == ".npz" else "npy",
            "tpm_hvg_obs_names": str(tpm_hvg_obs_path),
            "tpm_hvg_var_names": str(tpm_hvg_var_path),
        }
    if precomputed_tpm_stats is not None:
        _link_or_copy(precomputed_tpm_stats, tpm_stats_path)
        tpm_stats = _load_df_from_npz(tpm_stats_path)
        tpm_n_vars = int(tpm_stats.shape[0])
    else:
        _save_df_to_npz(tpm_stats, tpm_stats_path)
    hvgs_out.write_text("\n".join(hvgs), encoding="utf-8")
    nmf_params = _nmf_iter_params(k_values, n_iter, seed, iter_template)
    _save_df_to_npz(nmf_params, nmf_params_path)
    run_params = {
        "alpha_W": float(alpha_usage),
        "alpha_H": float(alpha_spectra),
        "l1_ratio": 0.0,
        "beta_loss": beta_loss,
        "solver": "cd" if beta_loss == "frobenius" else "mu",
        "tol": 1e-4,
        "max_iter": int(max_nmf_iter),
        "init": init,
    }
    with nmf_run_params_path.open("w", encoding="utf-8") as handle:
        yaml.dump(run_params, handle)

    manifest = {
        "engine": "fastcnmf.run_fast_prepare",
        "run_name": run_name,
        "output_dir": str(output_dir),
        "corrected_h5ad": str(corrected_h5ad),
        "tpm_h5ad": str(tpm_h5ad) if tpm_h5ad is not None else None,
        "hvgs_txt": str(hvgs_txt),
        "normalized_counts": str(normalized_counts_path) if h5ad_written else None,
        "normalized_counts_h5ad": str(normalized_counts_path) if h5ad_written else None,
        "normalized_counts_npy": str(normalized_counts_npy_path) if npy_written else None,
        "normalized_counts_obs_names": str(normalized_counts_obs_path) if npy_written else None,
        "normalized_counts_var_names": str(normalized_counts_var_path) if npy_written else None,
        "tpm": tpm_out,
        "tpm_link_mode": tpm_link_mode,
        "tpm_stats": str(tpm_stats_path),
        "nmf_params": str(nmf_params_path),
        "nmf_run_params": str(nmf_run_params_path),
        "hvgs_out": str(hvgs_out),
        "n_obs": int(norm_x.shape[0]),
        "n_vars": int(norm_x.shape[1]),
        "tpm_n_vars": int(tpm_n_vars),
        "k_values": list(map(int, k_values)),
        "n_iter": int(n_iter),
        "seed": int(seed),
        "norm_dtype": norm_dtype,
        "norm_store": norm_store,
        "precomputed_norm_npy": str(precomputed_norm_npy) if precomputed_norm_npy is not None else None,
        "precomputed_tpm_stats": str(precomputed_tpm_stats) if precomputed_tpm_stats is not None else None,
        **tpm_hvg_cache,
    }
    manifest_path = tmp_dir / f"{run_name}.fast_prepare_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
