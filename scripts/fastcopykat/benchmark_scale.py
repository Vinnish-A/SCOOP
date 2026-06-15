from __future__ import annotations

import argparse
import json
import resource
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd

from fastcopykat import FastCopyKatConfig, run_fastcopykat
from fastcopykat.io import write_copykat_outputs
from fastcopykat.reference import filter_copykat_chromosomes, normalize_bins, normalize_gene_annotation, read_table


def simulate_counts(
    annotation: pd.DataFrame,
    *,
    n_cells: int,
    n_genes: int,
    aneuploid_fraction: float,
    seed: int,
) -> tuple[pd.DataFrame, pd.Series]:
    rng = np.random.default_rng(seed)
    anno = annotation.drop_duplicates("gene").copy()
    per_chrom = max(25, n_genes // max(1, anno["chrom"].nunique()))
    selected = (
        anno.groupby("chrom", sort=False, group_keys=False)[["gene", "chrom", "start", "end", "chrompos", "abspos"]]
        .apply(lambda frame: frame.sample(min(len(frame), per_chrom), random_state=int(seed)))
        .sort_values(["chrom", "chrompos"])
    )
    if selected.shape[0] > n_genes:
        selected = selected.sample(n_genes, random_state=int(seed)).sort_values(["chrom", "chrompos"])
    genes = selected["gene"].to_numpy()
    chrom = selected["chrom"].astype(str).to_numpy()

    n_aneuploid = int(round(n_cells * aneuploid_fraction))
    n_diploid = n_cells - n_aneuploid
    labels = np.array(["diploid"] * n_diploid + ["aneuploid"] * n_aneuploid, dtype=object)
    cell_names = np.array([f"dip{i:05d}" for i in range(n_diploid)] + [f"tumor{i:05d}" for i in range(n_aneuploid)])

    gene_base = rng.gamma(shape=1.6, scale=2.5, size=genes.shape[0]) + 0.02
    library = rng.lognormal(mean=0.0, sigma=0.35, size=n_cells)
    lam = gene_base[:, None] * library[None, :]

    cnv = np.ones((genes.shape[0], n_cells), dtype=np.float32)
    tumor_slice = slice(n_diploid, n_cells)
    clone_switch = n_diploid + n_aneuploid // 2
    cnv[np.isin(chrom, ["7", "8", "20"]), tumor_slice] *= 1.9
    cnv[np.isin(chrom, ["10", "13", "14"]), tumor_slice] *= 0.45
    cnv[np.isin(chrom, ["1", "5"]), n_diploid:clone_switch] *= 1.55
    cnv[np.isin(chrom, ["3", "11"]), clone_switch:n_cells] *= 0.55
    lam *= cnv
    lam = np.clip(lam, 0, 150)

    counts = rng.poisson(lam).astype(np.uint16)
    counts_df = pd.DataFrame(counts, index=genes, columns=cell_names)
    truth = pd.Series(labels, index=cell_names, name="truth")
    return counts_df, truth


def score_predictions(prediction: pd.DataFrame, truth: pd.Series) -> dict:
    pred = prediction.set_index("cell.names")["copykat.pred"].reindex(truth.index)
    valid = pred.notna()
    truth = truth.loc[valid]
    pred = pred.loc[valid]
    labels = ("diploid", "aneuploid")
    per_label = {}
    recalls = []
    for label in labels:
        mask = truth == label
        recall = float((pred[mask] == label).mean()) if mask.any() else float("nan")
        per_label[label] = recall
        if not np.isnan(recall):
            recalls.append(recall)
    return {
        "accuracy": float((pred == truth).mean()),
        "balanced_accuracy": float(np.mean(recalls)),
        "per_label_recall": per_label,
        "truth_counts": truth.value_counts().to_dict(),
        "prediction_counts": pred.value_counts().to_dict(),
        "overlap_cells": int(valid.sum()),
    }


def run_one_size(
    *,
    annotation: pd.DataFrame,
    bins: pd.DataFrame | None,
    output_dir: Path,
    n_cells: int,
    n_genes: int,
    aneuploid_fraction: float,
    seed: int,
    write_outputs: bool,
    output_mode: str,
    save_inputs: bool,
) -> dict:
    counts, truth = simulate_counts(
        annotation,
        n_cells=n_cells,
        n_genes=n_genes,
        aneuploid_fraction=aneuploid_fraction,
        seed=seed,
    )
    cfg = FastCopyKatConfig(
        min_gene_per_cell=200,
        min_gene_per_chromosome=5,
        low_detection_rate=0.05,
        upper_detection_rate=0.10,
        window_size=25,
        min_cluster_cells=max(5, min(25, n_cells // 40)),
    )
    start = perf_counter()
    result = run_fastcopykat(
        counts,
        annotation,
        bins=bins,
        sample_name=f"synthetic_{n_cells}",
        config=cfg,
    )
    compute_seconds = perf_counter() - start
    write_seconds = 0.0
    outputs = {}
    if save_inputs:
        input_dir = output_dir / f"synthetic_{n_cells}"
        input_dir.mkdir(parents=True, exist_ok=True)
        counts.to_csv(input_dir / f"synthetic_{n_cells}_counts.tsv", sep="\t")
        truth.rename("truth").to_csv(input_dir / f"synthetic_{n_cells}_truth.tsv", sep="\t", header=True)
        outputs["counts"] = str(input_dir / f"synthetic_{n_cells}_counts.tsv")
        outputs["truth"] = str(input_dir / f"synthetic_{n_cells}_truth.tsv")
    if write_outputs:
        start = perf_counter()
        outputs.update(write_copykat_outputs(
            prediction=result.prediction,
            cna=result.cna,
            cell_scores=result.cell_scores,
            output_dir=output_dir / f"synthetic_{n_cells}",
            sample_name=f"synthetic_{n_cells}",
            manifest=result.manifest,
            mode=output_mode,
        ))
        write_seconds = perf_counter() - start
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return {
        "n_cells": n_cells,
        "n_genes_requested": n_genes,
        "n_genes_retained": result.manifest["n_retained_genes"],
        "n_bins": result.manifest["n_bins"],
        "compute_seconds": compute_seconds,
        "write_seconds": write_seconds,
        "total_seconds": compute_seconds + write_seconds,
        "max_rss_mb": usage.ru_maxrss / 1024.0,
        "stage_timings": result.manifest["timings"],
        "score": score_predictions(result.prediction, truth),
        "outputs": outputs,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run larger FastCopyKAT synthetic benchmarks.")
    parser.add_argument("--gene-annotation", required=True)
    parser.add_argument("--bins", default=None)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--sizes", type=int, nargs="+", default=[1000, 3000])
    parser.add_argument("--genes", type=int, default=12000)
    parser.add_argument("--aneuploid-fraction", type=float, default=0.70)
    parser.add_argument("--seed", type=int, default=20260615)
    parser.add_argument("--write-outputs", action="store_true", default=False)
    parser.add_argument(
        "--output-mode",
        choices=["compact", "parquet", "copykat-tsv"],
        default="compact",
        help="output mode used with --write-outputs",
    )
    parser.add_argument("--save-inputs", action="store_true", default=False)
    args = parser.parse_args()

    annotation = filter_copykat_chromosomes(normalize_gene_annotation(read_table(Path(args.gene_annotation))), genome="hg20")
    bins = None
    if args.bins:
        bins = filter_copykat_chromosomes(normalize_bins(read_table(Path(args.bins))), genome="hg20")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for offset, n_cells in enumerate(args.sizes):
        results.append(
            run_one_size(
                annotation=annotation,
                bins=bins,
                output_dir=output_dir,
                n_cells=n_cells,
                n_genes=args.genes,
                aneuploid_fraction=args.aneuploid_fraction,
                seed=args.seed + offset,
                write_outputs=args.write_outputs,
                output_mode=args.output_mode,
                save_inputs=args.save_inputs,
            )
        )
    payload = {
        "schema_version": "fastcopykat.synthetic_scale_benchmark.v1",
        "sizes": args.sizes,
        "n_genes": args.genes,
        "aneuploid_fraction": args.aneuploid_fraction,
        "write_outputs": args.write_outputs,
        "output_mode": args.output_mode,
        "save_inputs": args.save_inputs,
        "results": results,
    }
    out_json = output_dir / "fastcopykat_scale_benchmark.json"
    out_json.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
