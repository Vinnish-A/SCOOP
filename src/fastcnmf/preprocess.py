from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from .harmony2_compat import harmony2_moe_correct


def _save_df_to_npz(obj: pd.DataFrame, path: Path, *, compressed: bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = np.savez_compressed if compressed else np.savez
    writer(path, data=obj.values, index=obj.index.values, columns=obj.columns.values)


def _scale_zero_center_false_array(x: np.ndarray, *, output_dtype: str = "float32") -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    std = x.std(axis=0, ddof=1)
    std[std == 0] = 1.0
    x /= std
    if np.isnan(x).any():
        raise ValueError("NaNs produced while scaling normalized counts")
    if output_dtype == "float32":
        return x.astype(np.float32)
    if output_dtype == "float64":
        return x
    raise ValueError("output_dtype must be float32 or float64")


def _sparse_mean_var(matrix) -> tuple[np.ndarray, np.ndarray]:
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler(with_mean=False)
    scaler.fit(matrix)
    return np.asarray(scaler.mean_), np.asarray(scaler.var_)


def _dense_mean_var(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(matrix, dtype=np.float64)
    return x.mean(axis=0), x.var(axis=0, ddof=0)


def _write_matrix(matrix, path: Path) -> tuple[str, str]:
    path.parent.mkdir(parents=True, exist_ok=True)
    if sparse.issparse(matrix):
        out = Path(str(path) + ".npz")
        sparse.save_npz(out, matrix.tocsr())
        return str(out), "sparse_npz"
    out = Path(str(path) + ".npy")
    np.save(out, np.asarray(matrix))
    return str(out), "npy"


def stdscale_quantile_ceiling(adata, *, max_value: float | None = None, quantile_thresh: float | None = 0.9999) -> None:
    """Scale genes with cNMF-compatible zero-preserving quantile clipping."""

    import scanpy as sc
    from scipy import sparse

    sc.pp.scale(adata, zero_center=False, max_value=max_value)
    if quantile_thresh is None:
        return
    values = np.asarray(adata.X.todense() if sparse.issparse(adata.X) else adata.X)
    thresh = np.quantile(values.reshape(-1), quantile_thresh)
    adata.X[adata.X > thresh] = thresh


def _dense_x(adata) -> np.ndarray:
    from scipy import sparse

    return np.asarray(adata.X.todense() if sparse.issparse(adata.X) else adata.X, dtype=np.float64)


def run_fast_harmony_preprocess(
    *,
    input_h5ad: Path,
    output_prefix: Path,
    sample_key: str = "sample_id",
    n_top_genes: int = 3000,
    librarysize_targetsum: float = 1e4,
    theta: float = 1.0,
    lamb: float = 1.0,
    max_iter_harmony: int = 20,
    seed: int = 20260614,
    quantile_thresh: float = 0.9999,
    write_core_cache: bool = False,
    core_dtype: str = "float32",
    write_corrected_h5ad: bool = True,
    write_tp10k_h5ad: bool = True,
) -> dict:
    """Prepare cNMF-compatible inputs with FastCNMF-owned preprocessing."""

    import anndata as ad
    import scanpy as sc

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(input_h5ad)
    adata.var_names_make_unique()
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    if sample_key not in adata.obs.columns:
        raise ValueError(f"sample_key {sample_key!r} is not present in obs")

    tp10k = adata.copy()
    sc.pp.normalize_total(tp10k, target_sum=librarysize_targetsum, copy=False)

    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_top_genes)
    hvg_mask = adata.var["highly_variable"].to_numpy()
    if int(hvg_mask.sum()) == 0:
        raise ValueError("HVG selection produced zero genes")

    anorm = sc.pp.normalize_total(adata, target_sum=librarysize_targetsum, copy=True)
    anorm = anorm[:, hvg_mask].copy()
    stdscale_quantile_ceiling(anorm, quantile_thresh=quantile_thresh)
    sc.pp.pca(anorm, use_highly_variable=True, zero_center=True, random_state=seed)

    hvg = adata[:, hvg_mask].copy()
    stdscale_quantile_ceiling(hvg, quantile_thresh=quantile_thresh)

    result = harmony2_moe_correct(
        _dense_x(hvg),
        np.asarray(anorm.obsm["X_pca"], dtype=np.float64),
        hvg.obs,
        sample_key,
        lamb=lamb,
        theta=theta,
        max_iter_harmony=max_iter_harmony,
        ncores=1,
        verbose=False,
        random_state=seed,
    )

    corrected = ad.AnnData(
        X=result.x_corr,
        obs=hvg.obs.copy(),
        var=pd.DataFrame(index=hvg.var_names.copy()),
        obsm={
            "X_pca": np.asarray(anorm.obsm["X_pca"], dtype=np.float64),
            "X_pca_harmony": result.x_pca_harmony,
        },
    )
    corrected.layers["counts"] = corrected.X.copy()
    corrected.uns["fastcnmf_preprocess"] = {
        "method": "fastcnmf_harmony_moe",
        "sample_key": sample_key,
        "n_top_genes": int(n_top_genes),
        "librarysize_targetsum": float(librarysize_targetsum),
        "theta": float(theta),
        "lamb": float(lamb),
        "max_iter_harmony": int(max_iter_harmony),
        "seed": int(seed),
        "quantile_thresh": float(quantile_thresh),
    }

    corrected_path = Path(str(output_prefix) + ".Corrected.HVG.Varnorm.h5ad")
    tp10k_path = Path(str(output_prefix) + ".TP10K.h5ad")
    hvg_path = Path(str(output_prefix) + ".Corrected.HVGs.txt")
    norm_core_path = Path(str(output_prefix) + f".NormCounts.{core_dtype}.npy")
    norm_core_obs_path = Path(str(output_prefix) + ".NormCounts.obs.txt")
    norm_core_var_path = Path(str(output_prefix) + ".NormCounts.var.txt")
    tpm_stats_path = Path(str(output_prefix) + ".TP10K.stats.df.npz")
    tpm_hvg_raw_base = Path(str(output_prefix) + ".TP10K.HVG.raw")
    tpm_hvg_scaled_base = Path(str(output_prefix) + ".TP10K.HVG.scaled")
    tpm_hvg_obs_path = Path(str(output_prefix) + ".TP10K.HVG.obs.txt")
    tpm_hvg_var_path = Path(str(output_prefix) + ".TP10K.HVG.var.txt")
    if write_corrected_h5ad:
        corrected.write_h5ad(corrected_path)
    if write_tp10k_h5ad:
        tp10k.write_h5ad(tp10k_path)
    hvg_path.write_text("\n".join(map(str, hvg.var_names)) + "\n", encoding="utf-8")

    core_cache = {}
    if write_core_cache:
        norm_core = _scale_zero_center_false_array(np.asarray(result.x_corr, dtype=np.float64).copy(), output_dtype=core_dtype)
        np.save(norm_core_path, norm_core)
        norm_core_obs_path.write_text("\n".join(map(str, hvg.obs_names)), encoding="utf-8")
        norm_core_var_path.write_text("\n".join(map(str, hvg.var_names)), encoding="utf-8")

        if sparse.issparse(tp10k.X):
            mean, var = _sparse_mean_var(tp10k.X)
        else:
            mean, var = _dense_mean_var(tp10k.X)
        tpm_stats = pd.DataFrame(
            [mean, np.sqrt(var)],
            index=["__mean", "__std"],
            columns=tp10k.var.index,
        ).T
        _save_df_to_npz(tpm_stats, tpm_stats_path)

        tp10k_hvg = tp10k[:, hvg_mask].copy()
        raw_path, raw_format = _write_matrix(tp10k_hvg.X, tpm_hvg_raw_base)
        sc.pp.scale(tp10k_hvg, zero_center=False)
        scaled_path, scaled_format = _write_matrix(tp10k_hvg.X, tpm_hvg_scaled_base)
        tpm_hvg_obs_path.write_text("\n".join(map(str, tp10k_hvg.obs_names)), encoding="utf-8")
        tpm_hvg_var_path.write_text("\n".join(map(str, tp10k_hvg.var_names)), encoding="utf-8")
        core_cache = {
            "normalized_counts_npy": str(norm_core_path),
            "normalized_counts_obs_names": str(norm_core_obs_path),
            "normalized_counts_var_names": str(norm_core_var_path),
            "tpm_stats": str(tpm_stats_path),
            "tpm_hvg_raw": raw_path,
            "tpm_hvg_raw_format": raw_format,
            "tpm_hvg_scaled": scaled_path,
            "tpm_hvg_scaled_format": scaled_format,
            "tpm_hvg_obs_names": str(tpm_hvg_obs_path),
            "tpm_hvg_var_names": str(tpm_hvg_var_path),
            "core_dtype": core_dtype,
        }

    manifest = {
        "engine": "fastcnmf.run_fast_harmony_preprocess",
        "input_h5ad": str(input_h5ad),
        "output_prefix": str(output_prefix),
        "corrected_h5ad": str(corrected_path) if write_corrected_h5ad else None,
        "tp10k_h5ad": str(tp10k_path) if write_tp10k_h5ad else None,
        "hvg_txt": str(hvg_path),
        "sample_key": sample_key,
        "n_obs": int(corrected.n_obs),
        "n_vars": int(corrected.n_vars),
        "n_top_genes": int(n_top_genes),
        "librarysize_targetsum": float(librarysize_targetsum),
        "theta": float(theta),
        "lamb": float(lamb),
        "max_iter_harmony": int(max_iter_harmony),
        "seed": int(seed),
        "write_core_cache": bool(write_core_cache),
        "write_corrected_h5ad": bool(write_corrected_h5ad),
        "write_tp10k_h5ad": bool(write_tp10k_h5ad),
        **core_cache,
    }
    manifest_path = output_prefix.parent / "fastcnmf_preprocess_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
