from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


DEFAULT_DATASETS = (
    "public:h5ad/canonical/quick_test/public_O_GSE154795_10samples_1000cells.h5ad",
    "private:h5ad/canonical/quick_test/private_overall_sim_10samples_1000cells.h5ad",
)


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def parse_dataset(value: str) -> tuple[str, Path]:
    if ":" not in value:
        path = Path(value)
        return safe_name(path.stem), path
    dataset_id, path = value.split(":", 1)
    return safe_name(dataset_id), Path(path)


def export_sample_counts(
    *,
    h5ad_path: Path,
    dataset_id: str,
    sample_key: str,
    sample_id: str,
    n_cells: int,
    layer: str | None,
    seed: int,
    output_dir: Path,
) -> tuple[Path, dict]:
    obj = ad.read_h5ad(h5ad_path)
    sample_mask = obj.obs[sample_key].astype(str).to_numpy() == str(sample_id)
    sample_indices = np.flatnonzero(sample_mask)
    if sample_indices.size < n_cells:
        raise ValueError(f"sample {sample_id} has {sample_indices.size} cells, fewer than requested {n_cells}")
    rng = np.random.default_rng(seed)
    chosen = np.sort(rng.choice(sample_indices, size=n_cells, replace=False))
    sub = obj[chosen, :].copy()
    matrix = sub.layers[layer] if layer else (sub.layers["counts"] if "counts" in sub.layers else sub.X)
    if sp.issparse(matrix):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix)
    counts = pd.DataFrame(matrix.T, index=sub.var_names.astype(str), columns=sub.obs_names.astype(str))
    counts = counts.round().astype(np.uint32)

    sample_dir = output_dir / dataset_id / safe_name(sample_id)
    sample_dir.mkdir(parents=True, exist_ok=True)
    counts_path = sample_dir / "counts.tsv"
    counts.to_csv(counts_path, sep="\t")
    meta = {
        "dataset_id": dataset_id,
        "sample_id": sample_id,
        "h5ad_path": str(h5ad_path),
        "cells_available": int(sample_indices.size),
        "cells_exported": int(n_cells),
        "genes_exported": int(counts.shape[0]),
        "counts_path": str(counts_path),
    }
    (sample_dir / "sample_manifest.json").write_text(json.dumps(meta, indent=2, sort_keys=True), encoding="utf-8")
    return counts_path, meta


def parse_time_verbose(stderr: str) -> dict:
    metrics = {}
    patterns = {
        "user_seconds": r"User time \(seconds\): ([0-9.]+)",
        "system_seconds": r"System time \(seconds\): ([0-9.]+)",
        "elapsed": r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)",
        "max_rss_kb": r"Maximum resident set size \(kbytes\): ([0-9]+)",
        "file_outputs": r"File system outputs: ([0-9]+)",
        "exit_status": r"Exit status: ([0-9]+)",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, stderr)
        if not match:
            continue
        value = match.group(1)
        if key == "elapsed":
            metrics["wall_seconds"] = parse_elapsed(value)
        elif key == "max_rss_kb":
            metrics["max_rss_mb"] = int(value) / 1024.0
        elif key in {"file_outputs", "exit_status"}:
            metrics[key] = int(value)
        else:
            metrics[key] = float(value)
    return metrics


def parse_elapsed(value: str) -> float:
    parts = value.split(":")
    if len(parts) == 3:
        hours, minutes, seconds = parts
        return int(hours) * 3600 + int(minutes) * 60 + float(seconds)
    if len(parts) == 2:
        minutes, seconds = parts
        return int(minutes) * 60 + float(seconds)
    return float(value)


def run_timed(command: list[str], *, cwd: Path, env: dict[str, str] | None = None, timeout: int | None = None) -> dict:
    full_env = os.environ.copy()
    if env:
        full_env.update(env)
    proc = subprocess.run(
        ["/usr/bin/time", "-v", *command],
        cwd=cwd,
        env=full_env,
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    metrics = parse_time_verbose(proc.stderr)
    return {
        "command": command,
        "returncode": proc.returncode,
        "stdout": proc.stdout[-4000:],
        "stderr_tail": proc.stderr[-4000:],
        "metrics": metrics,
    }


def read_prediction(path: Path) -> pd.Series:
    frame = pd.read_csv(path, sep="\t")
    return frame.set_index("cell.names")["copykat.pred"]


def normalize_copykat_label(value: str) -> str:
    lowered = str(value).lower()
    if "not.defined" in lowered or lowered in {"nan", "none"}:
        return "not.defined"
    if "aneuploid" in lowered:
        return "aneuploid"
    if "diploid" in lowered:
        return "diploid"
    return str(value)


def compare_predictions(reference: Path, candidate: Path) -> dict:
    ref = read_prediction(reference)
    cand = read_prediction(candidate)
    merged = pd.concat([ref.rename("copykat_r"), cand.rename("fastcopykat")], axis=1, join="inner")
    confusion = pd.crosstab(merged["copykat_r"], merged["fastcopykat"])
    normalized = merged.apply(lambda col: col.map(normalize_copykat_label))
    defined = normalized["copykat_r"] != "not.defined"
    norm_confusion = pd.crosstab(normalized.loc[defined, "copykat_r"], normalized.loc[defined, "fastcopykat"])
    return {
        "overlap_cells": int(merged.shape[0]),
        "agreement": float((merged["copykat_r"] == merged["fastcopykat"]).mean()),
        "normalized_defined_cells": int(defined.sum()),
        "normalized_agreement_defined": float(
            (normalized.loc[defined, "copykat_r"] == normalized.loc[defined, "fastcopykat"]).mean()
        )
        if defined.any()
        else None,
        "copykat_counts": ref.value_counts().to_dict(),
        "fastcopykat_counts": cand.value_counts().to_dict(),
        "copykat_counts_normalized": normalized["copykat_r"].value_counts().to_dict(),
        "fastcopykat_counts_normalized": normalized["fastcopykat"].value_counts().to_dict(),
        "confusion": {
            str(row): {str(col): int(value) for col, value in confusion.loc[row].items()}
            for row in confusion.index
        },
        "confusion_normalized_defined": {
            str(row): {str(col): int(value) for col, value in norm_confusion.loc[row].items()}
            for row in norm_confusion.index
        },
    }


def output_size(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file())


def write_markdown(results: list[dict], output_md: Path) -> None:
    rows = []
    for item in results:
        fast = item["fastcopykat"]["metrics"]
        ref = item["copykat_r"]["metrics"]
        compare = item["comparison"]
        speedup = ref.get("wall_seconds", float("nan")) / fast.get("wall_seconds", float("nan"))
        rows.append(
            "| {dataset} | {sample} | {cells} | {fast_wall:.2f}s | {ref_wall:.2f}s | {speedup:.2f}x | "
            "{fast_rss:.2f} GB | {ref_rss:.2f} GB | {agreement:.2%} | {norm_agreement} | {fast_counts} | {ref_counts} |".format(
                dataset=item["dataset_id"],
                sample=item["sample_id"],
                cells=item["cells_exported"],
                fast_wall=fast.get("wall_seconds", float("nan")),
                ref_wall=ref.get("wall_seconds", float("nan")),
                speedup=speedup,
                fast_rss=fast.get("max_rss_mb", 0) / 1024.0,
                ref_rss=ref.get("max_rss_mb", 0) / 1024.0,
                agreement=compare["agreement"],
                norm_agreement=(
                    "{:.2%}".format(compare["normalized_agreement_defined"])
                    if compare.get("normalized_agreement_defined") is not None
                    else "n/a"
                ),
                fast_counts=json.dumps(compare["fastcopykat_counts"], sort_keys=True),
                ref_counts=json.dumps(compare["copykat_counts"], sort_keys=True),
            )
        )
    output_md.write_text(
        "\n".join(
            [
                "# FastCopyKAT Project Sample Benchmark",
                "",
                "Each row is one H5AD sample subset run independently through FastCopyKAT compact output and R CopyKAT.",
                "",
                "| Dataset | Sample | Cells | Fast wall | CopyKAT wall | Speedup | Fast RSS | CopyKAT RSS | Raw agreement | Defined normalized agreement | Fast counts | CopyKAT counts |",
                "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |",
                *rows,
                "",
                "Raw agreement compares exact strings. Defined normalized agreement maps CopyKAT low-confidence labels to diploid/aneuploid and excludes `not.defined`. These quick-test samples do not include a biological CNV truth label.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark FastCopyKAT and R CopyKAT per sample on project H5AD tests.")
    parser.add_argument("--dataset", action="append", default=None, help="dataset_id:path.h5ad")
    parser.add_argument("--gene-annotation", default="tmp/fastcopykat_fixture/copykat_full_anno_hg20.tsv")
    parser.add_argument("--bins", default="tmp/fastcopykat_fixture/copykat_DNA_hg20_bins.tsv")
    parser.add_argument("--output-dir", default="tmp/fastcopykat_project_sample_benchmark")
    parser.add_argument("--sample-key", default="sample_id")
    parser.add_argument("--samples-per-dataset", type=int, default=2)
    parser.add_argument("--cells-per-sample", type=int, default=300)
    parser.add_argument("--layer", default=None)
    parser.add_argument("--seed", type=int, default=20260615)
    parser.add_argument("--timeout", type=int, default=1800)
    args = parser.parse_args()

    cwd = Path.cwd()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    datasets = [parse_dataset(value) for value in (args.dataset or DEFAULT_DATASETS)]
    results = []
    for dataset_index, (dataset_id, h5ad_path) in enumerate(datasets):
        obj = ad.read_h5ad(h5ad_path, backed="r")
        counts = obj.obs[args.sample_key].astype(str).value_counts()
        selected_samples = list(counts[counts >= args.cells_per_sample].index[: args.samples_per_dataset])
        obj.file.close()
        for sample_index, sample_id in enumerate(selected_samples):
            sample_seed = args.seed + dataset_index * 100 + sample_index
            counts_path, sample_meta = export_sample_counts(
                h5ad_path=h5ad_path,
                dataset_id=dataset_id,
                sample_key=args.sample_key,
                sample_id=sample_id,
                n_cells=args.cells_per_sample,
                layer=args.layer,
                seed=sample_seed,
                output_dir=output_dir / "inputs",
            )
            run_root = output_dir / "runs" / dataset_id / safe_name(sample_id)
            fast_dir = run_root / "fastcopykat"
            copykat_dir = run_root / "copykat_r"
            fast = run_timed(
                [
                    sys.executable,
                    "-m",
                    "fastcopykat",
                    "run",
                    "--input",
                    str(counts_path),
                    "--gene-annotation",
                    args.gene_annotation,
                    "--bins",
                    args.bins,
                    "--output-dir",
                    str(fast_dir),
                    "--sample-name",
                    "sample",
                    "--output-mode",
                    "compact",
                    "--min-gene-per-cell",
                    "200",
                    "--min-gene-per-chromosome",
                    "5",
                ],
                cwd=cwd,
                env={"PYTHONPATH": "src"},
                timeout=args.timeout,
            )
            ref = run_timed(
                [
                    "Rscript",
                    "scripts/fastcopykat/run_copykat_reference_synthetic.R",
                    str(counts_path),
                    str(copykat_dir),
                    "sample_r",
                    "1",
                ],
                cwd=cwd,
                timeout=args.timeout,
            )
            fast_pred = fast_dir / "sample_copykat_prediction.txt"
            ref_pred = copykat_dir / "sample_r_copykat_prediction.txt"
            comparison = compare_predictions(ref_pred, fast_pred) if fast_pred.exists() and ref_pred.exists() else {}
            result = {
                **sample_meta,
                "fastcopykat": fast,
                "copykat_r": ref,
                "comparison": comparison,
                "output_sizes": {
                    "fastcopykat_bytes": output_size(fast_dir),
                    "copykat_r_bytes": output_size(copykat_dir),
                },
            }
            results.append(result)
            (run_root / "result.json").parent.mkdir(parents=True, exist_ok=True)
            (run_root / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")

    summary = {
        "schema_version": "fastcopykat.project_sample_benchmark.v1",
        "cells_per_sample": args.cells_per_sample,
        "samples_per_dataset": args.samples_per_dataset,
        "results": results,
    }
    summary_json = output_dir / "summary.json"
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(results, output_dir / "summary.md")
    print(json.dumps({"summary_json": str(summary_json), "summary_md": str(output_dir / "summary.md"), "runs": len(results)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
