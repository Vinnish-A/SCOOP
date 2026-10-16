"""Checkpointed process-parallel sweeps with OpenMP inside NetworKit Leiden."""

from __future__ import annotations

import hashlib
import json
import multiprocessing
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy import sparse

_GRAPH = None


def _initialize(graph_path: str, threads: int) -> None:
    import networkit as nk

    global _GRAPH
    from threadpoolctl import threadpool_limits

    threadpool_limits(limits=1, user_api="blas")
    nk.setNumberOfThreads(threads)
    matrix = sparse.load_npz(graph_path)
    upper = sparse.triu(matrix, k=1, format="coo")
    _GRAPH = nk.GraphFromCoo(upper, n=matrix.shape[0], weighted=True, directed=False)


def _cluster(resolution: float, seed: int, iterations: int) -> tuple[np.ndarray, dict]:
    import networkit as nk

    nk.setSeed(seed, True)
    started = time.perf_counter()
    cpu_started = time.process_time()
    algorithm = nk.community.ParallelLeiden(_GRAPH, randomize=True, iterations=iterations, gamma=resolution)
    algorithm.run()
    labels = np.asarray(algorithm.getPartition().getVector(), dtype=np.int64)
    _, labels = np.unique(labels, return_inverse=True)
    elapsed = time.perf_counter() - started
    cpu_seconds = time.process_time() - cpu_started
    return labels, {
        "resolution": resolution,
        "seed": seed,
        "iterations": iterations,
        "threads": nk.getMaxNumberOfThreads(),
        "pid": os.getpid(),
        "elapsed_seconds": elapsed,
        "cpu_seconds": cpu_seconds,
        "cpu_to_wall_ratio": cpu_seconds / max(elapsed, 1e-9),
        "n_clusters": int(np.unique(labels).size),
        "backend": "networkit_parallel_leiden",
        "networkit_version": nk.__version__,
        "objective": "generalized_modularity",
        "iteration_policy": "fixed_budget",
    }


def run_jobs(adjacency, jobs, config):
    """Yield (resolution, seed, labels) as each independently seeded job completes."""
    import networkit as nk

    if not hasattr(nk.community, "ParallelLeiden"):
        raise RuntimeError("Managed runtime lacks NetworKit ParallelLeiden; repair it before clustering.")
    threads = int(config.get("threads_per_worker", 1))
    iterations = int(config.get("iterations", 2))
    if min(threads, iterations) < 1:
        raise ValueError("threads_per_worker and iterations must be positive")
    requested_workers = config.get("workers", "auto")
    if requested_workers != "auto" and int(requested_workers) < 1:
        raise ValueError("workers must be positive or 'auto'")
    root = Path(config["checkpoint_dir"])
    root.mkdir(parents=True, exist_ok=True)
    matrix = sparse.csr_matrix(adjacency, dtype=np.float64)
    matrix.sort_indices()
    if matrix.shape[0] != matrix.shape[1] or not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
        raise ValueError("ParallelLeiden requires a square, finite, nonnegative adjacency matrix")
    difference = matrix - matrix.T
    if difference.nnz and np.max(np.abs(difference.data)) > 1e-7:
        raise ValueError("ParallelLeiden requires an undirected symmetric graph")
    digest = hashlib.sha256()
    for value in (np.asarray(matrix.shape), matrix.indptr, matrix.indices, matrix.data):
        digest.update(value.tobytes())
    graph_hash = digest.hexdigest()
    graph_path = root / f"graph_{graph_hash}.npz"
    if not graph_path.exists():
        temporary = graph_path.with_suffix(".tmp.npz")
        sparse.save_npz(temporary, matrix, compressed=False)
        temporary.replace(graph_path)
    identity = dict(graph_sha256=graph_hash, iterations=iterations, threads=threads, networkit_version=nk.__version__)
    completed = 0
    pending = []
    for resolution, seed in jobs:
        path = root / f"r{resolution:g}_s{seed}.npz"
        metadata_path = path.with_suffix(".json")
        if path.exists() and metadata_path.exists():
            record = json.loads(metadata_path.read_text())
            if all(record.get(k) == v for k, v in identity.items()):
                labels = np.load(path)["labels"]
                if len(labels) != matrix.shape[0]:
                    raise ValueError(f"Checkpoint cell count mismatch: {path}")
                completed += 1
                print(f"LEIDEN_RESUMED {completed}/{len(jobs)} resolution={resolution} seed={seed}", flush=True)
                yield resolution, seed, labels
                continue
            raise ValueError(f"Checkpoint settings changed; use a new checkpoint directory: {path}")
        pending.append((resolution, seed))
    if not pending:
        return
    # Bound task concurrency by the current graph, memory, CPU affinity and load.
    cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    load = os.getloadavg()[0] if hasattr(os, "getloadavg") else 0.0
    free_cpus = max(1, int(cpus - load))
    available_memory = os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    memory_per_worker = 6 * (matrix.data.nbytes + matrix.indices.nbytes + matrix.indptr.nbytes) + 256 * 1024**2
    worker_limit = 4 if requested_workers == "auto" else int(requested_workers)
    workers = max(1, min(worker_limit, len(pending), free_cpus // threads, available_memory // memory_per_worker))
    print(
        f"LEIDEN_PARALLEL workers={workers} threads_per_worker={threads} jobs={len(pending)} nodes={matrix.shape[0]}",
        flush=True,
    )
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=multiprocessing.get_context("spawn"),
        initializer=_initialize,
        initargs=(str(graph_path), threads),
    ) as pool:
        futures = {pool.submit(_cluster, float(res), int(seed), iterations): (res, seed) for res, seed in pending}
        for future in as_completed(futures):
            resolution, seed = futures[future]
            labels, record = future.result()
            record.update(identity)
            record["workers"] = workers
            record["requested_workers"] = requested_workers
            record["available_cpus_after_load"] = free_cpus
            record["available_memory_bytes"] = available_memory
            record["completed_jobs"] = completed + 1
            record["total_jobs"] = len(jobs)
            path = root / f"r{resolution:g}_s{seed}.npz"
            temporary = path.with_suffix(".tmp.npz")
            np.savez_compressed(temporary, labels=labels)
            temporary.replace(path)
            metadata_path = path.with_suffix(".json")
            temporary_meta = metadata_path.with_suffix(".tmp.json")
            temporary_meta.write_text(json.dumps(record, indent=2) + "\n")
            temporary_meta.replace(metadata_path)
            completed += 1
            print(
                f"LEIDEN_COMPLETE {completed}/{len(jobs)} resolution={resolution} seed={seed} "
                f"clusters={record['n_clusters']} seconds={record['elapsed_seconds']:.2f} "
                f"cpu/wall={record['cpu_to_wall_ratio']:.2f}",
                flush=True,
            )
            yield resolution, seed, labels
