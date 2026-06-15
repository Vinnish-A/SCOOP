from __future__ import annotations

import multiprocessing as mp
import os
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


_X = None
_VAR_NAMES: np.ndarray | None = None
_RUN_PARAMS: dict | None = None
_OUTPUT_TEMPLATE: str | None = None
_COMPRESSED = False
_BACKEND = "exact"
_MINIBATCH_BATCH_SIZE = 1024
_MU_DTYPE = "float32"
_MU_EPS = 1e-8


@dataclass(frozen=True)
class FactorizeSummary:
    run_name: str
    output_dir: str
    jobs_total: int
    jobs_completed: int
    workers: int
    scheduler: str


def _load_df_from_npz(path: Path) -> pd.DataFrame:
    with np.load(path, allow_pickle=True) as f:
        return pd.DataFrame(**f)


def _save_df_to_npz(obj: pd.DataFrame, path: Path, *, compressed: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = np.savez_compressed if compressed else np.savez
    writer(path, data=obj.values, index=obj.index.values, columns=obj.columns.values)


def _cnmf_tmp(output_dir: Path, run_name: str) -> Path:
    return output_dir / run_name / "cnmf_tmp"


def _fast_prepare_manifest_path(output_dir: Path, run_name: str) -> Path:
    return _cnmf_tmp(output_dir, run_name) / f"{run_name}.fast_prepare_manifest.json"


def _read_names(path: Path) -> np.ndarray:
    return np.array(path.read_text(encoding="utf-8").splitlines(), dtype=object)


def _iter_spectra_template(output_dir: Path, run_name: str) -> str:
    return str(_cnmf_tmp(output_dir, run_name) / f"{run_name}.spectra.k_%d.iter_%d.df.npz")


def _load_nmf_kwargs(output_dir: Path, run_name: str) -> dict:
    import yaml

    path = _cnmf_tmp(output_dir, run_name) / f"{run_name}.nmf_idvrun_params.yaml"
    with path.open("r", encoding="utf-8") as handle:
        return yaml.load(handle, Loader=yaml.FullLoader)


def _load_replicate_params(output_dir: Path, run_name: str) -> pd.DataFrame:
    path = _cnmf_tmp(output_dir, run_name) / f"{run_name}.nmf_params.df.npz"
    return _load_df_from_npz(path)


def _load_matrix(output_dir: Path, run_name: str):
    manifest_path = _fast_prepare_manifest_path(output_dir, run_name)
    if manifest_path.is_file():
        import json

        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        npy = manifest.get("normalized_counts_npy")
        var_names = manifest.get("normalized_counts_var_names")
        if npy and var_names and Path(npy).is_file() and Path(var_names).is_file():
            return np.load(npy, mmap_mode="r"), _read_names(Path(var_names))

    path = _cnmf_tmp(output_dir, run_name) / f"{run_name}.norm_counts.h5ad"
    try:
        import scanpy as sc

        adata = sc.read_h5ad(path)
    except ModuleNotFoundError:
        import anndata as ad

        adata = ad.read_h5ad(path)
    return adata.X, adata.var.index.to_numpy()


def _pending_jobs(params: pd.DataFrame, output_template: str, skip_completed: bool) -> list[dict]:
    jobs: list[dict] = []
    for idx, row in params.iterrows():
        k = int(row["n_components"])
        iteration = int(row["iter"])
        output = Path(output_template % (k, iteration))
        if skip_completed and output.exists():
            continue
        jobs.append(
            {
                "index": int(idx),
                "k": k,
                "iteration": iteration,
                "seed": int(row["nmf_seed"]),
            }
        )
    return sorted(jobs, key=lambda job: (-job["k"], job["iteration"], job["index"]))


def _run_one_job(job: dict) -> dict:
    if _X is None or _VAR_NAMES is None or _RUN_PARAMS is None or _OUTPUT_TEMPLATE is None:
        raise RuntimeError("factorize worker was not initialized")
    started = time.perf_counter()
    kwargs = dict(_RUN_PARAMS)
    if _BACKEND == "exact":
        from sklearn.decomposition import non_negative_factorization

        kwargs["random_state"] = job["seed"]
        kwargs["n_components"] = job["k"]
        usages, spectra, n_iter = non_negative_factorization(_X, **kwargs)
        del usages
    elif _BACKEND == "minibatch":
        from sklearn.decomposition import MiniBatchNMF

        model = MiniBatchNMF(
            n_components=job["k"],
            init=kwargs.get("init", "random"),
            batch_size=_MINIBATCH_BATCH_SIZE,
            beta_loss=kwargs.get("beta_loss", "frobenius"),
            tol=kwargs.get("tol", 1e-4),
            max_iter=kwargs.get("max_iter", 200),
            alpha_W=kwargs.get("alpha_W", 0.0),
            alpha_H=kwargs.get("alpha_H", 0.0),
            l1_ratio=kwargs.get("l1_ratio", 0.0),
            random_state=job["seed"],
        )
        model.fit_transform(_X)
        spectra = model.components_
        n_iter = getattr(model, "n_iter_", kwargs.get("max_iter", 0))
    elif _BACKEND == "cupy-mu":
        spectra, n_iter = _run_cupy_mu_job(job, kwargs)
    else:
        raise ValueError(f"unsupported factorize backend: {_BACKEND}")
    spectra_df = pd.DataFrame(
        spectra,
        index=np.arange(1, job["k"] + 1),
        columns=_VAR_NAMES,
    )
    output = Path(_OUTPUT_TEMPLATE % (job["k"], job["iteration"]))
    _save_df_to_npz(spectra_df, output, compressed=_COMPRESSED)
    elapsed = time.perf_counter() - started
    return {
        "index": job["index"],
        "k": job["k"],
        "iteration": job["iteration"],
        "seed": job["seed"],
        "n_iter": int(n_iter),
        "elapsed_seconds": float(elapsed),
        "backend": _BACKEND,
        "output": str(output),
        "pid": os.getpid(),
    }


def _run_cupy_mu_job(job: dict, kwargs: dict) -> tuple[np.ndarray, int]:
    if kwargs.get("beta_loss", "frobenius") != "frobenius":
        raise ValueError("cupy-mu backend currently supports only frobenius beta_loss")
    if kwargs.get("alpha_W", 0.0) != 0.0 or kwargs.get("alpha_H", 0.0) != 0.0:
        raise ValueError("cupy-mu backend currently supports only unregularized NMF")

    import cupy as cp
    from sklearn.decomposition._nmf import _initialize_nmf

    x_cpu = np.asarray(_X, dtype=np.float32 if _MU_DTYPE == "float32" else np.float64)
    if x_cpu.min() < 0:
        raise ValueError("cupy-mu backend requires non-negative input")
    w_cpu, h_cpu = _initialize_nmf(
        x_cpu,
        n_components=job["k"],
        init=kwargs.get("init", "random"),
        random_state=job["seed"],
    )
    dtype = cp.float32 if _MU_DTYPE == "float32" else cp.float64
    x = cp.asarray(x_cpu, dtype=dtype)
    w = cp.asarray(w_cpu, dtype=dtype)
    h = cp.asarray(h_cpu, dtype=dtype)
    eps = cp.asarray(_MU_EPS, dtype=dtype)
    max_iter = int(kwargs.get("max_iter", 200))

    for _ in range(max_iter):
        h *= (w.T @ x) / cp.maximum((w.T @ w) @ h, eps)
        w *= (x @ h.T) / cp.maximum(w @ (h @ h.T), eps)

    cp.cuda.Stream.null.synchronize()
    spectra = cp.asnumpy(h).astype(np.float64, copy=False)
    del x, w, h
    cp.get_default_memory_pool().free_all_blocks()
    return spectra, max_iter


def run_fast_factorize(
    *,
    output_dir: Path,
    run_name: str,
    workers: int,
    skip_completed: bool = False,
    compressed: bool = False,
    backend: str = "exact",
    minibatch_batch_size: int = 1024,
    mu_dtype: str = "float32",
) -> dict:
    """Run cNMF-compatible NMF replicates with FastCNMF-owned scheduling.

    The implementation keeps cNMF's prepared files and output schema, but
    replaces cNMF's static worker-index sharding with a dynamic multiprocessing
    scheduler. On Linux the normalized matrix is loaded once in the parent and
    inherited by forked workers, avoiding repeated disk reads and reducing
    cross-worker memory pressure for read-only input.
    """

    if workers < 1:
        raise ValueError("workers must be positive")
    if backend not in {"exact", "minibatch", "cupy-mu"}:
        raise ValueError(f"unsupported backend: {backend}")
    if minibatch_batch_size < 1:
        raise ValueError("minibatch_batch_size must be positive")
    if mu_dtype not in {"float32", "float64"}:
        raise ValueError("mu_dtype must be float32 or float64")
    if backend == "cupy-mu" and workers != 1:
        workers = 1

    global _X, _VAR_NAMES, _RUN_PARAMS, _OUTPUT_TEMPLATE, _COMPRESSED, _BACKEND, _MINIBATCH_BATCH_SIZE, _MU_DTYPE
    _X, _VAR_NAMES = _load_matrix(output_dir, run_name)
    _RUN_PARAMS = _load_nmf_kwargs(output_dir, run_name)
    _OUTPUT_TEMPLATE = _iter_spectra_template(output_dir, run_name)
    _COMPRESSED = compressed
    _BACKEND = backend
    _MINIBATCH_BATCH_SIZE = minibatch_batch_size
    _MU_DTYPE = mu_dtype

    params = _load_replicate_params(output_dir, run_name)
    jobs = _pending_jobs(params, _OUTPUT_TEMPLATE, skip_completed=skip_completed)
    results: list[dict] = []
    if not jobs:
        summary = FactorizeSummary(
            run_name=run_name,
            output_dir=str(output_dir),
            jobs_total=int(len(params)),
            jobs_completed=0,
            workers=workers,
            scheduler=f"dynamic-fork-largest-k-first:{backend}",
        )
        return {**summary.__dict__, "backend": backend, "results": results}

    if workers == 1:
        results = [_run_one_job(job) for job in jobs]
    else:
        ctx = mp.get_context("fork")
        with ctx.Pool(processes=workers) as pool:
            for result in pool.imap_unordered(_run_one_job, jobs, chunksize=1):
                results.append(result)

    summary = FactorizeSummary(
        run_name=run_name,
        output_dir=str(output_dir),
        jobs_total=int(len(params)),
        jobs_completed=len(results),
        workers=workers,
        scheduler=f"dynamic-fork-largest-k-first:{backend}",
    )
    return {
        **summary.__dict__,
        "backend": backend,
        "minibatch_batch_size": minibatch_batch_size if backend == "minibatch" else None,
        "mu_dtype": mu_dtype if backend == "cupy-mu" else None,
        "results": sorted(results, key=lambda item: item["index"]),
    }
