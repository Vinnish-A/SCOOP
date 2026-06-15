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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from fastde.deseq2 import run_deseq2_wald
from fastde.io import read_pseudobulk_dir, write_json


def make_fixture(outdir: Path, *, n_genes: int, n_per_group: int, n_de: int, seed: int) -> Path:
    rng = np.random.default_rng(seed)
    ct_dir = outdir / "pseudobulk" / "T_cell"
    ct_dir.mkdir(parents=True, exist_ok=True)
    genes = [f"G{i:05d}" for i in range(n_genes)]
    base_mu = rng.gamma(shape=2.0, scale=80.0, size=n_genes)
    effect = np.zeros(n_genes)
    effect[:n_de] = np.log(2.0)
    rows = []
    meta = []
    for condition in ["ctrl", "test"]:
        for i in range(n_per_group):
            lib = rng.lognormal(mean=0.0, sigma=0.25)
            mu = base_mu * lib * np.exp(effect if condition == "test" else 0.0)
            alpha = 0.12
            lam = rng.gamma(shape=1.0 / alpha, scale=mu * alpha)
            rows.append(rng.poisson(lam))
            pid = f"T_cell|{condition}{i}|{condition}"
            meta.append({"pseudobulk_id": pid, "sample_id": f"{condition}{i}", "condition": condition, "cell_type": "T_cell", "n_cells": 100})
    counts = pd.DataFrame(rows, columns=genes)
    counts.insert(0, "pseudobulk_id", [m["pseudobulk_id"] for m in meta])
    counts.to_csv(ct_dir / "counts.tsv", sep="\t", index=False)
    pd.DataFrame(meta).to_csv(ct_dir / "metadata.tsv", sep="\t", index=False)
    return ct_dir


def run_fastde(indir: Path) -> tuple[pd.DataFrame, float, dict]:
    counts, meta = read_pseudobulk_dir(indir)
    start = perf_counter()
    result = run_deseq2_wald(counts, meta, condition_col="condition", ctrl_group="ctrl", test_group="test")
    seconds = perf_counter() - start
    outdir = indir.parent.parent / "contrasts" / "test_vs_ctrl" / indir.name
    outdir.mkdir(parents=True, exist_ok=True)
    result.table.to_csv(outdir / "de_fastde_deseq2.tsv", sep="\t", index=False)
    write_json(result.manifest, outdir / "fastde_deseq2_manifest.json")
    return result.table, seconds, result.manifest


def run_r_deseq2(indir: Path) -> tuple[pd.DataFrame | None, float | None, str | None]:
    check = subprocess.run(["Rscript", "-e", "cat(requireNamespace('DESeq2', quietly=TRUE))"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if check.stdout.strip() != "TRUE":
        return None, None, "DESeq2 is not installed"
    cmd = ["Rscript", str(ROOT / "r" / "run_pseudobulk_deseq2.R"), str(indir), "condition", "ctrl", "test"]
    start = perf_counter()
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    seconds = perf_counter() - start
    if proc.returncode != 0:
        return None, seconds, proc.stderr.strip()[-1000:]
    table = pd.read_csv(indir.parent.parent / "contrasts" / "test_vs_ctrl" / indir.name / "de_DESeq2.tsv", sep="\t")
    return table, seconds, None


def compare(candidate: pd.DataFrame, reference: pd.DataFrame, *, top_n: int) -> dict:
    cand = candidate.rename(columns={"log2FoldChange": "log2FoldChange"}).set_index("gene")
    ref = reference.set_index("gene")
    common = cand.index.intersection(ref.index)
    ref_p = "pvalue"
    cand_p = "pvalue"
    return {
        "n_common_genes": int(len(common)),
        "spearman_log2fc": float(cand.loc[common, "log2FoldChange"].corr(ref.loc[common, "log2FoldChange"], method="spearman")),
        "sign_agreement": float((np.sign(cand.loc[common, "log2FoldChange"]) == np.sign(ref.loc[common, "log2FoldChange"])).mean()),
        "top_n": int(top_n),
        "top_n_overlap_fraction": float(len(set(cand.sort_values(cand_p).head(top_n).index) & set(ref.sort_values(ref_p).head(top_n).index)) / top_n),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", default="tmp/fastde_deseq2_benchmark")
    ap.add_argument("--n-genes", type=int, default=10000)
    ap.add_argument("--n-per-group", type=int, default=6)
    ap.add_argument("--n-de", type=int, default=300)
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    indir = make_fixture(outdir, n_genes=args.n_genes, n_per_group=args.n_per_group, n_de=args.n_de, seed=args.seed)
    fast_table, fast_seconds, manifest = run_fastde(indir)
    ref_table, ref_seconds, ref_error = run_r_deseq2(indir)
    report = {
        "schema_version": "fastde.deseq2_benchmark.v1",
        "fixture": {"n_genes": args.n_genes, "n_per_group": args.n_per_group, "n_de": args.n_de, "seed": args.seed},
        "fastde_seconds": round(fast_seconds, 6),
        "fastde_manifest": manifest,
    }
    if ref_table is None:
        report["r_deseq2"] = {"status": "skipped_or_failed", "error": ref_error, "seconds": ref_seconds}
    else:
        report["r_deseq2_seconds"] = round(ref_seconds, 6)
        report["speedup_vs_r_deseq2"] = round(ref_seconds / fast_seconds, 6) if fast_seconds > 0 else None
        report["consistency"] = compare(fast_table, ref_table, top_n=min(100, args.n_genes))
    write_json(report, outdir / "benchmark_report.json")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
