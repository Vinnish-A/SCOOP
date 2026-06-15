from __future__ import annotations

import argparse
import json
from pathlib import Path

from .config import FastCopyKatConfig
from .core import run_fastcopykat
from .io import read_counts, write_copykat_outputs
from .reference import read_table


def cmd_run(args: argparse.Namespace) -> int:
    counts = read_counts(Path(args.input), layer=args.layer, transpose=args.transpose)
    annotation = read_table(Path(args.gene_annotation))
    bins = read_table(Path(args.bins)) if args.bins else None
    normals = _read_names(Path(args.normal_cells)) if args.normal_cells else []
    config = FastCopyKatConfig(
        genome=args.genome,
        min_gene_per_cell=args.min_gene_per_cell,
        min_gene_per_chromosome=args.min_gene_per_chromosome,
        low_detection_rate=args.low_detection_rate,
        upper_detection_rate=args.upper_detection_rate,
        window_size=args.window_size,
        bin_size=args.bin_size,
        segmentation_threshold=args.segmentation_threshold,
        min_cluster_cells=args.min_cluster_cells,
        max_baseline_clusters=args.max_baseline_clusters,
        distance=args.distance,
    )
    result = run_fastcopykat(
        counts,
        annotation,
        bins=bins,
        normal_cell_names=normals,
        sample_name=args.sample_name,
        config=config,
    )
    outputs = write_copykat_outputs(
        prediction=result.prediction,
        cna=result.cna,
        cell_scores=result.cell_scores,
        output_dir=Path(args.output_dir),
        sample_name=args.sample_name,
        manifest=result.manifest,
        mode=args.output_mode,
    )
    print(json.dumps(outputs, indent=2))
    return 0


def _read_names(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="fastcopykat")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="run Python CopyKAT-compatible CNV inference")
    run.add_argument("--input", required=True, help="gene-by-cell counts table or h5ad")
    run.add_argument("--gene-annotation", required=True, help="gene coordinate table")
    run.add_argument("--bins", default=None, help="optional CopyKAT-style bin coordinate table")
    run.add_argument("--output-dir", required=True)
    run.add_argument("--sample-name", default="sample")
    run.add_argument("--layer", default=None, help="h5ad layer to use")
    run.add_argument("--transpose", action="store_true", default=False, help="transpose non-h5ad input after reading")
    run.add_argument("--normal-cells", default=None, help="optional newline-delimited normal cell names")
    run.add_argument("--genome", default="hg20")
    run.add_argument("--min-gene-per-cell", type=int, default=200)
    run.add_argument("--min-gene-per-chromosome", type=int, default=5)
    run.add_argument("--low-detection-rate", type=float, default=0.05)
    run.add_argument("--upper-detection-rate", type=float, default=0.10)
    run.add_argument("--window-size", type=int, default=25)
    run.add_argument("--bin-size", type=int, default=220_000)
    run.add_argument("--segmentation-threshold", type=float, default=0.10)
    run.add_argument("--min-cluster-cells", type=int, default=5)
    run.add_argument("--max-baseline-clusters", type=int, default=6)
    run.add_argument("--distance", choices=["euclidean", "correlation"], default="euclidean")
    run.add_argument(
        "--output-mode",
        choices=["compact", "parquet", "copykat-tsv"],
        default="compact",
        help="compact writes prediction/scores/manifest only; parquet writes CNA as parquet; copykat-tsv writes the wide TSV",
    )
    run.set_defaults(func=cmd_run)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)
