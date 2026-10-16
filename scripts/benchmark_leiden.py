#!/usr/bin/env python
"""Measure real Leiden sweep throughput and quality on existing expression graphs."""

# ruff: noqa: E402 -- numerical imports must follow worker thread limits
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from time import perf_counter

# Independent workers should spend their CPU budget on Leiden rather than BLAS
# pools spawned by NumPy/sklearn imports. Set before any numerical imports.
for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import adjusted_rand_score

from scsp_agent_sop.parallel_leiden import run_jobs

PROFILES = {
    "reference_w1_t1_i10": (1, 1, 10),
    "w1_t1_i2": (1, 1, 2),
    "w4_t1_i2": (4, 1, 2),
    "w4_t2_i2": (4, 2, 2),
    "w4_t4_i2": (4, 4, 2),
    "w8_t1_i2": (8, 1, 2),
    "w2_t4_i2": (2, 4, 2),
    "w4_t2_i5": (4, 2, 5),
    "w4_t2_i10": (4, 2, 10),
}


def read_graph(path: Path):
    with h5py.File(path, "r") as file:
        group = file["obsp/connectivities_identity"]
        shape = tuple(group.attrs["shape"])
        constructor = sparse.csc_matrix if group.attrs["encoding-type"] == "csc_matrix" else sparse.csr_matrix
        return constructor((group["data"][:], group["indices"][:], group["indptr"][:]), shape=shape).tocsr()


def quality(matrix, labels, resolution):
    coo = matrix.tocoo()
    degree = np.asarray(matrix.sum(axis=1)).ravel()
    total = degree.sum()
    within = coo.data[labels[coo.row] == labels[coo.col]].sum() / total
    volumes = np.bincount(labels, weights=degree) / total
    return float(within - resolution * np.square(volumes).sum())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--profiles", nargs="+", choices=list(PROFILES), default=list(PROFILES))
    args = parser.parse_args()
    if (args.output / "summary.tsv").exists():
        raise FileExistsError("Choose a fresh benchmark output directory; cached jobs cannot measure execution speed.")
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    jobs = [(r, s) for r in (0.25, 0.75, 1.25) for s in (0, 1)]
    metadata = {
        "inputs": list(map(str, args.input)),
        "jobs": jobs,
        "logical_cpus": os.cpu_count(),
        "load_average": os.getloadavg(),
        "objective": "weighted undirected generalized modularity (gamma=resolution)",
        "iterations_note": (
            "NetworKit fixed iteration budget; differs from leidenalg convergence and igraph n_iterations."
        ),
        "quality_thresholds": {"median_ari_vs_reference": 0.75, "median_objective_loss": 0.01},
        "timing_includes": "pool startup, worker imports, graph loading, computation, checkpoints",
    }
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    # Charge the one-off main-process import separately. Every timed profile
    # still includes its own fresh worker/pool startup.
    import_started = perf_counter()
    import networkit

    metadata["networkit_version"] = networkit.__version__
    metadata["main_networkit_import_seconds"] = perf_counter() - import_started
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    for path in args.input:
        graph = read_graph(path)
        reference = {}
        reference_q = {}
        baseline_seconds = None
        graph_dir = args.output / path.parent.parent.name
        for name in ["reference_w1_t1_i10", *[p for p in args.profiles if p != "reference_w1_t1_i10"]]:
            workers, threads, iterations = PROFILES[name]
            config = {
                "workers": workers,
                "threads_per_worker": threads,
                "iterations": iterations,
                "checkpoint_dir": str(graph_dir / name),
            }
            started = perf_counter()
            labels = {(r, s): y for r, s, y in run_jobs(graph, jobs, config)}
            elapsed = perf_counter() - started
            if baseline_seconds is None:
                baseline_seconds = elapsed
                reference = labels
                reference_q = {j: quality(graph, y, j[0]) for j, y in labels.items()}
            for job, y in labels.items():
                q = quality(graph, y, job[0])
                rows.append(
                    {
                        "graph": path.parent.parent.name,
                        "nodes": graph.shape[0],
                        "nnz": graph.nnz,
                        "profile": name,
                        "workers": workers,
                        "threads_per_worker": threads,
                        "iterations": iterations,
                        "resolution": job[0],
                        "seed": job[1],
                        "n_clusters": int(np.unique(y).size),
                        "sweep_seconds": elapsed,
                        "speedup": baseline_seconds / elapsed,
                        "ari_vs_reference": adjusted_rand_score(reference[job], y),
                        "objective": q,
                        "objective_loss": reference_q[job] - q,
                        "seed_stability_ari": adjusted_rand_score(labels[(job[0], 0)], labels[(job[0], 1)]),
                    }
                )
            frame = pd.DataFrame(rows)
            frame.to_csv(args.output / "jobs.tsv", sep="\t", index=False)
            summary = frame.groupby(["graph", "profile"], sort=False).agg(
                nodes=("nodes", "first"),
                workers=("workers", "first"),
                threads=("threads_per_worker", "first"),
                iterations=("iterations", "first"),
                seconds=("sweep_seconds", "first"),
                speedup=("speedup", "first"),
                median_ari=("ari_vs_reference", "median"),
                min_ari=("ari_vs_reference", "min"),
                median_objective_loss=("objective_loss", "median"),
                median_seed_ari=("seed_stability_ari", "median"),
            )
            summary.to_csv(args.output / "summary.tsv", sep="\t")
            print("BENCHMARK", path.parent.parent.name, name, f"seconds={elapsed:.3f}", flush=True)


if __name__ == "__main__":
    main()
