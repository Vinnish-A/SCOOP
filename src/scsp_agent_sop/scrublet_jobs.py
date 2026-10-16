"""Independent library Scrublet jobs with count-bound resumable results."""
from __future__ import annotations

import hashlib
import json
import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from importlib.metadata import version
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

import numpy as np
from scipy import sparse


def _fit(path, parameters, threads):
    from threadpoolctl import threadpool_limits
    import anndata as ad
    import numba
    import scanpy as sc

    numba.set_num_threads(threads)
    start = perf_counter()
    with threadpool_limits(limits=threads):
        obj = ad.AnnData(sparse.load_npz(path))
        sc.pp.scrublet(obj, **parameters, threshold=None)
    return (
        obj.obs.doublet_score.to_numpy(float),
        obj.obs.predicted_doublet.to_numpy(bool),
        perf_counter() - start,
    )


def scrublet_jobs(adata, sample_key, counts_layer, parameters, *, workers="auto", threads=1, checkpoint_dir=None):
    if threads < 1 or (workers != "auto" and int(workers) < 1):
        raise ValueError("Scrublet workers and threads must be positive")
    if adata.obs[sample_key].isna().any():
        raise ValueError("Scrublet requires a library ID for every cell")
    if not adata.obs_names.is_unique:
        raise ValueError("Scrublet requires unique cell IDs")
    with TemporaryDirectory(prefix="scoop-scrublet-") as temporary:
        root = Path(checkpoint_dir) if checkpoint_dir else Path(temporary) / "results"
        root.mkdir(parents=True, exist_ok=True)
        pending, completed, largest = [], 0, 1
        groups = adata.obs.groupby(sample_key, observed=True, sort=False).indices
        total = len(groups)
        for sample, positions in groups.items():
            matrix = sparse.csr_matrix(adata.layers[counts_layer][positions]).copy()
            matrix.sum_duplicates()
            matrix.sort_indices()
            identity = {
                "sample": str(sample), "parameters": parameters,
                "versions": {name: version(name) for name in ("scanpy", "numpy", "scipy")},
                "shape": matrix.shape,
            }
            digest = hashlib.sha256(json.dumps(identity, sort_keys=True).encode())
            for array in (matrix.indptr, matrix.indices, matrix.data):
                digest.update(array.dtype.str.encode())
                digest.update(memoryview(array).cast("B"))
            digest.update(json.dumps(adata.obs_names[positions].tolist()).encode())
            digest.update(json.dumps(adata.var_names.tolist()).encode())
            key = digest.hexdigest()
            result_path = root / f"{key}.npz"
            if result_path.exists():
                with np.load(result_path, allow_pickle=False) as saved:
                    scores, pred = saved["scores"], saved["pred"]
                    record = json.loads(str(saved["record"]))
                if scores.shape != (len(positions),) or pred.shape != scores.shape:
                    raise ValueError(f"invalid Scrublet checkpoint: {result_path}")
                completed += 1
                print(f"SCRUBLET_RESUMED {completed}/{total} sample={sample}", flush=True)
                yield sample, positions, scores, pred, dict(record, resumed=True)
                continue
            input_path = Path(temporary) / f"{key}.npz"
            sparse.save_npz(input_path, matrix)
            largest = max(largest, matrix.data.nbytes + matrix.indices.nbytes + matrix.indptr.nbytes)
            pending.append((sample, positions, input_path, result_path))
        if not pending:
            return
        cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
        load = os.getloadavg()[0] if hasattr(os, "getloadavg") else 0
        free_cpus = max(1, int(cpus - load))
        memory = os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
        limit = 4 if workers == "auto" else int(workers)
        active = max(1, min(limit, len(pending), free_cpus // threads, memory // (largest * 12 + 1024**3)))
        print(f"SCRUBLET_PARALLEL workers={active} threads={threads} pending={len(pending)}/{total}", flush=True)
        with ProcessPoolExecutor(max_workers=active, mp_context=multiprocessing.get_context("spawn")) as pool:
            futures = {pool.submit(_fit, item[2], parameters, threads): item for item in pending}
            for future in as_completed(futures):
                sample, positions, input_path, result_path = futures[future]
                scores, pred, elapsed = future.result()
                if not np.isfinite(scores).all():
                    raise ValueError(f"Scrublet produced non-finite scores for library {sample}")
                completed += 1
                record = {
                    "backend": "scanpy_scrublet", "workers": active, "threads_per_worker": threads,
                    "elapsed_seconds": elapsed, "completed_jobs": completed, "total_jobs": total,
                    "resumed": False, "sample": str(sample), "parameters": parameters,
                }
                temp_result = result_path.with_suffix(".tmp.npz")
                np.savez_compressed(temp_result, scores=scores, pred=pred, record=json.dumps(record))
                temp_result.replace(result_path)
                input_path.unlink()
                print(f"SCRUBLET_DONE {completed}/{total} sample={sample} elapsed={elapsed:.2f}s", flush=True)
                yield sample, positions, scores, pred, record
