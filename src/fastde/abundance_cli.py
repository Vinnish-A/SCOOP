from __future__ import annotations

import argparse
import json

from .abundance import run_abundance


def cmd_abundance(args: argparse.Namespace) -> int:
    result = run_abundance(
        mode=args.mode,
        input_h5ad=args.input_h5ad,
        counts=args.counts,
        metadata=args.metadata,
        sample_key=args.sample_key,
        celltype_key=args.celltype_key,
        label_col=args.label_col,
        positive_label=args.positive_label,
        negative_label=args.negative_label,
        reference_level=args.reference_level,
        time_col=args.time_col,
        event_col=args.event_col,
        value_col=args.value_col,
        covariates=args.covariates,
        min_cells_per_sample=args.min_cells_per_sample,
        min_total_cells_per_celltype=args.min_total_cells_per_celltype,
        learning_rate=args.learning_rate,
        weight_decay=args.weight_decay,
        max_epochs=args.max_epochs,
        random_seed=args.random_seed,
        hidden_dim=args.hidden_dim,
        dropout=args.dropout,
        max_instances_per_sample=args.max_instances_per_sample,
        survival_loss=args.survival_loss,
        output_dir=args.output_dir,
    )
    print(json.dumps(result.manifest["outputs"], indent=2))
    return 0


def add_abundance_subparser(subparsers: argparse._SubParsersAction) -> None:
    abundance = subparsers.add_parser("abundance", help="run sample-level differential abundance and outcome association")
    abundance.add_argument("--mode", choices=["survival", "binary", "multiclass", "continuous", "condition"], required=True)
    abundance.add_argument("--input-h5ad", default=None)
    abundance.add_argument("--counts", default=None)
    abundance.add_argument("--metadata", default=None)
    abundance.add_argument("--sample-key", default="sample_id")
    abundance.add_argument("--celltype-key", default="cell_type_lvl3")
    abundance.add_argument("--label-col", default=None)
    abundance.add_argument("--positive-label", default=None)
    abundance.add_argument("--negative-label", default=None)
    abundance.add_argument("--reference-level", default=None)
    abundance.add_argument("--time-col", default=None)
    abundance.add_argument("--event-col", default=None)
    abundance.add_argument("--value-col", default=None)
    abundance.add_argument("--covariates", default=None)
    abundance.add_argument("--output-dir", required=True)
    abundance.add_argument("--min-cells-per-sample", type=int, default=20)
    abundance.add_argument("--min-total-cells-per-celltype", type=int, default=50)
    abundance.add_argument("--learning-rate", type=float, default=1e-3)
    abundance.add_argument("--weight-decay", type=float, default=1e-4)
    abundance.add_argument("--max-epochs", type=int, default=500)
    abundance.add_argument("--random-seed", type=int, default=0)
    abundance.add_argument("--hidden-dim", type=int, default=64)
    abundance.add_argument("--dropout", type=float, default=0.1)
    abundance.add_argument("--max-instances-per-sample", type=int, default=2000)
    abundance.add_argument("--survival-loss", choices=["cox", "cox_rank", "cox_plus_rank"], default="cox")
    abundance.set_defaults(func=cmd_abundance)
