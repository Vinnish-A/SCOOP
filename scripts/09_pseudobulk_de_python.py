#!/usr/bin/env python
"""Run Python-native pseudobulk condition DE for one cell-type directory."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fastde.deseq2 import run_deseq2_wald
from fastde.io import read_pseudobulk_dir
from scsp_agent_sop.storage import ensure_dir, write_json


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("pseudobulk_dir")
    ap.add_argument("condition_col")
    ap.add_argument("ctrl_group")
    ap.add_argument("test_group")
    ap.add_argument("--method", choices=["deseq2-wald"], default="deseq2-wald")
    ap.add_argument("--min-total-count", type=int, default=15)
    ap.add_argument("--min-samples-per-group", type=int, default=2)
    ap.add_argument("--max-iter", type=int, default=50)
    ap.add_argument("--beta-tol", type=float, default=1e-8)
    args = ap.parse_args()

    indir = Path(args.pseudobulk_dir)
    counts, metadata = read_pseudobulk_dir(indir)
    result = run_deseq2_wald(
        counts,
        metadata,
        condition_col=args.condition_col,
        ctrl_group=args.ctrl_group,
        test_group=args.test_group,
        min_total_count=args.min_total_count,
        min_samples_per_group=args.min_samples_per_group,
        max_iter=args.max_iter,
        beta_tol=args.beta_tol,
    )
    outdir = indir.parent.parent / "contrasts" / f"{args.test_group}_vs_{args.ctrl_group}" / indir.name
    ensure_dir(outdir)
    result.table.to_csv(outdir / "de_fastde_deseq2.tsv", sep="\t", index=False)
    result.design.to_csv(outdir / "design_matrix_fastde.tsv", sep="\t")
    result.size_factors.to_csv(outdir / "size_factors_fastde.tsv", sep="\t", header=True)
    result.dispersions.to_csv(outdir / "dispersions_fastde.tsv", sep="\t", header=True)
    write_json(result.manifest, outdir / "fastde_deseq2_manifest.json")


if __name__ == "__main__":
    main()
