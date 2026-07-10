#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from scsp_agent_sop.de import compare_de_tables, read_pseudobulk_dir, run_pseudobulk_logcpm_welch
from scsp_agent_sop.storage import ensure_dir, write_json


def make_fixture(outdir: Path, *, n_genes: int, n_per_group: int, n_de: int, seed: int) -> Path:
    rng = np.random.default_rng(seed)
    ct_dir = outdir / "pseudobulk" / "T_cell"
    ensure_dir(ct_dir)
    genes = [f"G{i:05d}" for i in range(n_genes)]
    samples = []
    meta = []
    base_mu = rng.gamma(shape=2.0, scale=80.0, size=n_genes)
    de = np.zeros(n_genes)
    de[:n_de] = np.log(2.0)
    for group in ["ctrl", "test"]:
        for i in range(n_per_group):
            sample = f"T_cell|{group}{i}|{group}"
            lib = rng.lognormal(mean=0.0, sigma=0.25)
            mu = base_mu * lib * np.exp(de if group == "test" else 0.0)
            alpha = 0.12
            shape = 1.0 / alpha
            scale = mu * alpha
            lam = rng.gamma(shape=shape, scale=scale)
            samples.append(rng.poisson(lam))
            meta.append({"pseudobulk_id": sample, "cell_type": "T_cell", "sample_id": f"{group}{i}", "condition": group, "n_cells": 100})
    counts = pd.DataFrame(samples, columns=genes)
    counts.insert(0, "pseudobulk_id", [m["pseudobulk_id"] for m in meta])
    counts.to_csv(ct_dir / "counts.tsv", sep="\t", index=False)
    pd.DataFrame(meta).to_csv(ct_dir / "metadata.tsv", sep="\t", index=False)
    return ct_dir


def run_python(indir: Path) -> tuple[pd.DataFrame, float, dict]:
    counts, meta = read_pseudobulk_dir(indir)
    start = perf_counter()
    result = run_pseudobulk_logcpm_welch(counts, meta, condition_col="condition", ctrl_group="ctrl", test_group="test")
    seconds = perf_counter() - start
    outdir = indir.parent.parent / "contrasts" / "test_vs_ctrl" / indir.name
    ensure_dir(outdir)
    result.table.to_csv(outdir / "de_python_logcpm_welch.tsv", sep="\t", index=False)
    write_json(result.manifest, outdir / "de_python_manifest.json")
    return result.table, seconds, result.manifest


def run_edger(indir: Path) -> tuple[pd.DataFrame, float]:
    cmd = ["Rscript", str(ROOT / "r" / "run_pseudobulk_edger.R"), str(indir), "condition", "ctrl", "test"]
    start = perf_counter()
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    seconds = perf_counter() - start
    out = indir.parent.parent / "contrasts" / "test_vs_ctrl" / indir.name / "de_edgeR.tsv"
    return pd.read_csv(out, sep="\t"), seconds


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", default="tmp/pseudobulk_de_python_benchmark")
    ap.add_argument("--n-genes", type=int, default=3000)
    ap.add_argument("--n-per-group", type=int, default=4)
    ap.add_argument("--n-de", type=int, default=120)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--skip-edger", action="store_true")
    args = ap.parse_args()

    outdir = Path(args.output_dir)
    ensure_dir(outdir)
    indir = make_fixture(outdir, n_genes=args.n_genes, n_per_group=args.n_per_group, n_de=args.n_de, seed=args.seed)
    py_table, py_seconds, py_manifest = run_python(indir)
    report = {
        "schema_version": "scoop.pseudobulk_de_python_benchmark.v1",
        "fixture": {
            "n_genes": args.n_genes,
            "n_per_group": args.n_per_group,
            "n_de": args.n_de,
            "seed": args.seed,
        },
        "python_seconds": round(py_seconds, 6),
        "python_manifest": py_manifest,
    }
    if not args.skip_edger:
        edge_table, edge_seconds = run_edger(indir)
        report["edgeR_seconds"] = round(edge_seconds, 6)
        report["speedup_vs_edgeR"] = round(edge_seconds / py_seconds, 6) if py_seconds > 0 else None
        report["consistency"] = compare_de_tables(py_table, edge_table, top_n=min(100, args.n_genes))
    write_json(report, outdir / "benchmark_report.json")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
