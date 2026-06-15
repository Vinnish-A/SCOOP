from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad

from .config import FastCNVConfig
from .core import run_fastcnv, run_fastcnv_anndata, run_fastcnv_pooled_anndata
from .io import read_counts, read_obs, write_outputs, write_pooled_outputs
from .reference import read_gene_metadata


def cmd_run(args: argparse.Namespace) -> int:
    input_path = Path(args.input)
    gene_metadata = read_gene_metadata(Path(args.gene_metadata))
    reference_label = args.reference_label
    if reference_label and "," in reference_label:
        reference_label = [item for item in reference_label.split(",") if item]
    cfg = FastCNVConfig(
        scale_on_reference_label=not args.no_scale_on_reference_label,
        threshold_percentile=args.threshold_percentile,
        window_size=args.window_size,
        window_step=args.window_step,
        top_n_genes=args.top_n_genes,
        cluster_k=args.cluster_k,
        cluster_h=args.cluster_h,
        merge_cnv=not args.no_merge_cnv,
        merge_threshold=args.merge_threshold,
    )
    if input_path.suffix.lower() == ".h5ad" and not args.transpose:
        adata = ad.read_h5ad(input_path)
        obs = read_obs(Path(args.obs) if args.obs else None)
        result = run_fastcnv_anndata(
            adata,
            gene_metadata,
            layer=args.layer,
            obs=obs,
            reference_var=args.reference_var,
            reference_label=reference_label,
            config=cfg,
            sample_name=args.sample_name,
            compute_clusters=not args.no_clusters,
            compute_classification=not args.no_classification,
            densify_all=args.h5ad_mode == "dense",
        )
    else:
        counts = read_counts(input_path, layer=args.layer, transpose=args.transpose)
        obs = read_obs(Path(args.obs) if args.obs else None, h5ad_path=input_path if input_path.suffix.lower() == ".h5ad" else None)
        result = run_fastcnv(
            counts,
            gene_metadata,
            obs=obs,
            reference_var=args.reference_var,
            reference_label=reference_label,
            config=cfg,
            sample_name=args.sample_name,
            compute_clusters=not args.no_clusters,
            compute_classification=not args.no_classification,
        )
    outputs = write_outputs(result=result, output_dir=Path(args.output_dir), sample_name=args.sample_name, mode=args.output_mode)
    print(json.dumps(outputs, indent=2))
    return 0


def cmd_run_pooled(args: argparse.Namespace) -> int:
    input_path = Path(args.input)
    if input_path.suffix.lower() != ".h5ad":
        raise ValueError("run-pooled currently expects a merged H5AD input")
    gene_metadata = read_gene_metadata(Path(args.gene_metadata))
    reference_label = args.reference_label
    if reference_label and "," in reference_label:
        reference_label = [item for item in reference_label.split(",") if item]
    cfg = FastCNVConfig(
        scale_on_reference_label=not args.no_scale_on_reference_label,
        threshold_percentile=args.threshold_percentile,
        window_size=args.window_size,
        window_step=args.window_step,
        top_n_genes=args.top_n_genes,
        cluster_k=args.cluster_k,
        cluster_h=args.cluster_h,
        merge_cnv=not args.no_merge_cnv,
        merge_threshold=args.merge_threshold,
    )
    adata = ad.read_h5ad(input_path)
    obs = read_obs(Path(args.obs) if args.obs else None)
    result = run_fastcnv_pooled_anndata(
        adata,
        gene_metadata,
        sample_key=args.sample_key,
        layer=args.layer,
        obs=obs,
        reference_var=args.reference_var,
        reference_label=reference_label,
        config=cfg,
        sample_name=args.sample_name,
        compute_clusters=not args.no_clusters,
        compute_classification=not args.no_classification,
        densify_all=args.h5ad_mode == "dense",
        n_jobs=args.n_jobs,
        min_reference_cells_per_sample=args.min_reference_cells_per_sample,
    )
    outputs = write_pooled_outputs(result=result, output_dir=Path(args.output_dir), sample_name=args.sample_name, mode=args.output_mode)
    print(json.dumps(outputs, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="fastcnvpy")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="run a Python fastCNV-compatible CNVCalling pipeline")
    run.add_argument("--input", required=True, help="gene-by-cell table or h5ad")
    run.add_argument("--gene-metadata", required=True, help="fastCNV geneMetadata TSV/CSV/Parquet")
    run.add_argument("--output-dir", required=True)
    run.add_argument("--sample-name", default="sample")
    run.add_argument("--obs", default=None, help="optional cell metadata table")
    run.add_argument("--reference-var", default=None)
    run.add_argument("--reference-label", default=None, help="single label or comma-separated labels")
    run.add_argument("--layer", default=None)
    run.add_argument("--transpose", action="store_true", default=False)
    run.add_argument("--h5ad-mode", choices=["dense", "sparse"], default="dense", help="dense is fastest; sparse lowers peak memory")
    run.add_argument("--threshold-percentile", type=float, default=0.01)
    run.add_argument("--window-size", type=int, default=150)
    run.add_argument("--window-step", type=int, default=10)
    run.add_argument("--top-n-genes", type=int, default=7000)
    run.add_argument("--no-scale-on-reference-label", action="store_true", default=False)
    run.add_argument("--cluster-k", type=int, default=None)
    run.add_argument("--cluster-h", type=float, default=None)
    run.add_argument("--no-clusters", action="store_true", default=False)
    run.add_argument("--no-merge-cnv", action="store_true", default=False)
    run.add_argument("--merge-threshold", type=float, default=0.98)
    run.add_argument("--no-classification", action="store_true", default=False)
    run.add_argument("--output-mode", choices=["compact", "parquet", "tsv"], default="compact")
    run.set_defaults(func=cmd_run)

    pooled = sub.add_parser("run-pooled", help="run fastCNV pooled-reference workflow on a merged H5AD split by sample")
    pooled.add_argument("--input", required=True, help="merged h5ad")
    pooled.add_argument("--gene-metadata", required=True, help="fastCNV geneMetadata TSV/CSV/Parquet")
    pooled.add_argument("--output-dir", required=True)
    pooled.add_argument("--sample-name", default="pooled")
    pooled.add_argument("--sample-key", default="sample_id")
    pooled.add_argument("--obs", default=None, help="optional cell metadata table")
    pooled.add_argument("--reference-var", default=None)
    pooled.add_argument("--reference-label", default=None, help="single label or comma-separated labels")
    pooled.add_argument("--min-reference-cells-per-sample", type=int, default=5)
    pooled.add_argument("--layer", default=None)
    pooled.add_argument("--h5ad-mode", choices=["dense", "sparse"], default="dense", help="dense is fastest; sparse lowers peak memory")
    pooled.add_argument("--n-jobs", type=int, default=1, help="threaded per-sample score workers; keep low to cap memory")
    pooled.add_argument("--threshold-percentile", type=float, default=0.01)
    pooled.add_argument("--window-size", type=int, default=150)
    pooled.add_argument("--window-step", type=int, default=10)
    pooled.add_argument("--top-n-genes", type=int, default=7000)
    pooled.add_argument("--no-scale-on-reference-label", action="store_true", default=False)
    pooled.add_argument("--cluster-k", type=int, default=None)
    pooled.add_argument("--cluster-h", type=float, default=None)
    pooled.add_argument("--no-clusters", action="store_true", default=False)
    pooled.add_argument("--no-merge-cnv", action="store_true", default=False)
    pooled.add_argument("--merge-threshold", type=float, default=0.98)
    pooled.add_argument("--no-classification", action="store_true", default=False)
    pooled.add_argument("--output-mode", choices=["compact", "parquet", "tsv"], default="compact")
    pooled.set_defaults(func=cmd_run_pooled)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)
