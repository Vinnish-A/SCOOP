#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import resource
import sys
from pathlib import Path
from time import perf_counter

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from fastde.io import write_json
from fastde.markers import run_cosg_markers_dense, run_cosg_markers_sparse


def _rss_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def _select_matrix(adata: ad.AnnData, layer: str | None):
    return adata.layers[layer] if layer else adata.X


def _subset_adata(adata: ad.AnnData, *, max_cells: int | None, max_genes: int | None, seed: int) -> ad.AnnData:
    if (max_cells is None or adata.n_obs <= max_cells) and (max_genes is None or adata.n_vars <= max_genes):
        return adata
    rng = np.random.default_rng(seed)
    obs_idx = np.arange(adata.n_obs)
    var_idx = np.arange(adata.n_vars)
    if max_cells is not None and adata.n_obs > max_cells:
        obs_idx = np.sort(rng.choice(obs_idx, size=max_cells, replace=False))
    if max_genes is not None and adata.n_vars > max_genes:
        var_idx = np.sort(rng.choice(var_idx, size=max_genes, replace=False))
    return adata[obs_idx, var_idx].copy()


def _run_worker(args: dict) -> dict:
    start = perf_counter()
    before = _rss_mb()
    adata = ad.read_h5ad(args["input"])
    if args["groupby"] not in adata.obs:
        raise KeyError(f"groupby {args['groupby']!r} is not present in obs")
    adata = _subset_adata(adata, max_cells=args["max_cells"], max_genes=args["max_genes"], seed=args["seed"])
    matrix = _select_matrix(adata, args["layer"])
    labels = adata.obs[args["groupby"]]
    if args["backend"] == "sparse":
        result = run_cosg_markers_sparse(matrix, labels, adata.var_names, top_n=args["top_n"], min_pct=args["min_pct"])
    elif args["backend"] == "dense":
        result = run_cosg_markers_dense(matrix, labels, adata.var_names, top_n=args["top_n"], min_pct=args["min_pct"])
    else:
        raise ValueError(args["backend"])
    outdir = Path(args["output_dir"])
    outdir.mkdir(parents=True, exist_ok=True)
    table_path = outdir / f"markers_{args['backend']}.tsv"
    result.table.to_csv(table_path, sep="\t", index=False)
    manifest_path = write_json(result.manifest, outdir / f"markers_{args['backend']}_manifest.json")
    elapsed = perf_counter() - start
    return {
        "backend": args["backend"],
        "seconds_total": round(elapsed, 6),
        "seconds_algorithm": result.manifest["seconds"],
        "rss_before_mb": round(before, 3),
        "rss_peak_mb": round(_rss_mb(), 3),
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "n_groups": int(adata.obs[args["groupby"]].nunique()),
        "matrix_type": type(matrix).__name__,
        "table": str(table_path),
        "manifest": str(manifest_path),
    }


def _run_isolated(args: dict) -> dict:
    ctx = mp.get_context("spawn")
    with ctx.Pool(1) as pool:
        return pool.apply(_run_worker, (args,))


def _compare_tables(left_path: str, right_path: str, *, top_n: int) -> dict:
    left = pd.read_csv(left_path, sep="\t")
    right = pd.read_csv(right_path, sep="\t")
    groups = sorted(set(left["group"].astype(str)) & set(right["group"].astype(str)))
    overlaps = []
    exact_rank = []
    for group in groups:
        l = left[left["group"].astype(str) == group].sort_values("rank").head(top_n)
        r = right[right["group"].astype(str) == group].sort_values("rank").head(top_n)
        overlaps.append(len(set(l["gene"]) & set(r["gene"])) / max(1, min(top_n, len(l), len(r))))
        common_rank = l[["rank", "gene"]].merge(r[["rank", "gene"]], on=["rank", "gene"])
        exact_rank.append(len(common_rank) / max(1, min(top_n, len(l), len(r))))
    return {
        "n_groups_compared": int(len(groups)),
        "mean_top_overlap": float(np.mean(overlaps)) if overlaps else None,
        "min_top_overlap": float(np.min(overlaps)) if overlaps else None,
        "mean_exact_rank_match": float(np.mean(exact_rank)) if exact_rank else None,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="h5ad/canonical/quick_test/public_O_GSE154795_24samples_3000cells_balanced.h5ad")
    ap.add_argument("--groupby", default="sample_id")
    ap.add_argument("--layer", default="counts")
    ap.add_argument("--output-dir", default="tmp/fastde_cosg_benchmark")
    ap.add_argument("--top-n", type=int, default=50)
    ap.add_argument("--min-pct", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--dense-max-cells", type=int, default=6000)
    ap.add_argument("--dense-max-genes", type=int, default=5000)
    ap.add_argument("--full-sparse", action="store_true")
    ap.add_argument("--only-full-sparse", action="store_true")
    args = ap.parse_args()

    base = {
        "input": args.input,
        "groupby": args.groupby,
        "layer": args.layer,
        "output_dir": args.output_dir,
        "top_n": args.top_n,
        "min_pct": args.min_pct,
        "seed": args.seed,
    }
    dense_args = {**base, "backend": "dense", "max_cells": args.dense_max_cells, "max_genes": args.dense_max_genes, "output_dir": str(Path(args.output_dir) / "subset_dense")}
    sparse_subset_args = {**base, "backend": "sparse", "max_cells": args.dense_max_cells, "max_genes": args.dense_max_genes, "output_dir": str(Path(args.output_dir) / "subset_sparse")}
    sparse_full_args = {**base, "backend": "sparse", "max_cells": None, "max_genes": None, "output_dir": str(Path(args.output_dir) / "full_sparse")}

    report = {
        "schema_version": "fastde.cosg_benchmark.v1",
        "input": os.path.abspath(args.input),
        "groupby": args.groupby,
        "layer": args.layer,
        "top_n": args.top_n,
        "min_pct": args.min_pct,
        "subset": {"max_cells": args.dense_max_cells, "max_genes": args.dense_max_genes, "seed": args.seed},
    }
    if not args.only_full_sparse:
        report["subset_dense"] = _run_isolated(dense_args)
        report["subset_sparse"] = _run_isolated(sparse_subset_args)
        report["subset_consistency"] = _compare_tables(report["subset_dense"]["table"], report["subset_sparse"]["table"], top_n=args.top_n)
    if args.full_sparse or args.only_full_sparse:
        report["full_sparse"] = _run_isolated(sparse_full_args)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    write_json(report, out / "benchmark_report.json")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
