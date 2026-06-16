from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad

from .deseq2 import run_deseq2_wald
from .io import read_pseudobulk_dir, write_json
from .markers import run_cosg_markers, run_wilcoxon_markers
from .abundance_cli import add_abundance_subparser


def cmd_deseq2(args: argparse.Namespace) -> int:
    counts, meta = read_pseudobulk_dir(args.pseudobulk_dir)
    result = run_deseq2_wald(
        counts,
        meta,
        condition_col=args.condition_col,
        ctrl_group=args.ctrl_group,
        test_group=args.test_group,
        min_total_count=args.min_total_count,
        min_samples_per_group=args.min_samples_per_group,
        max_iter=args.max_iter,
        beta_tol=args.beta_tol,
    )
    outdir = Path(args.output_dir) if args.output_dir else Path(args.pseudobulk_dir).parent.parent / "contrasts" / f"{args.test_group}_vs_{args.ctrl_group}" / Path(args.pseudobulk_dir).name
    outdir.mkdir(parents=True, exist_ok=True)
    table_path = outdir / "de_fastde_deseq2.tsv"
    result.table.to_csv(table_path, sep="\t", index=False)
    result.design.to_csv(outdir / "design_matrix_fastde.tsv", sep="\t")
    result.size_factors.to_csv(outdir / "size_factors_fastde.tsv", sep="\t", header=True)
    result.dispersions.to_csv(outdir / "dispersions_fastde.tsv", sep="\t", header=True)
    manifest_path = write_json(result.manifest, outdir / "fastde_deseq2_manifest.json")
    print(json.dumps({"table": str(table_path), "manifest": str(manifest_path)}, indent=2))
    return 0


def cmd_markers(args: argparse.Namespace) -> int:
    adata = ad.read_h5ad(args.input)
    if args.groupby not in adata.obs.columns:
        raise KeyError(f"groupby {args.groupby!r} is not present in obs")
    matrix = adata.layers[args.layer] if args.layer else adata.X
    if args.method == "cosg":
        result = run_cosg_markers(matrix, adata.obs[args.groupby], adata.var_names, top_n=args.top_n, min_pct=args.min_pct)
    else:
        result = run_wilcoxon_markers(matrix, adata.obs[args.groupby], adata.var_names, top_n=args.top_n, min_pct=args.min_pct)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    table_path = outdir / f"markers_{args.method}.tsv"
    result.table.to_csv(table_path, sep="\t", index=False)
    manifest_path = write_json(result.manifest, outdir / f"markers_{args.method}_manifest.json")
    print(json.dumps({"table": str(table_path), "manifest": str(manifest_path)}, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="fastde")
    sub = parser.add_subparsers(dest="command", required=True)

    de = sub.add_parser("deseq2", help="run Python DESeq2-like pseudobulk NB Wald DE")
    de.add_argument("pseudobulk_dir")
    de.add_argument("condition_col")
    de.add_argument("ctrl_group")
    de.add_argument("test_group")
    de.add_argument("--output-dir", default=None)
    de.add_argument("--min-total-count", type=int, default=10)
    de.add_argument("--min-samples-per-group", type=int, default=2)
    de.add_argument("--max-iter", type=int, default=50)
    de.add_argument("--beta-tol", type=float, default=1e-8)
    de.set_defaults(func=cmd_deseq2)

    markers = sub.add_parser("markers", help="run marker-gene detection")
    markers.add_argument("--input", required=True)
    markers.add_argument("--groupby", required=True)
    markers.add_argument("--output-dir", required=True)
    markers.add_argument("--layer", default=None)
    markers.add_argument("--method", choices=["cosg", "wilcoxon"], default="cosg")
    markers.add_argument("--top-n", type=int, default=100)
    markers.add_argument("--min-pct", type=float, default=0.05)
    markers.set_defaults(func=cmd_markers)

    add_abundance_subparser(sub)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)
