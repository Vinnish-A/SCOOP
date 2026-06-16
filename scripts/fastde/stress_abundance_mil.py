from __future__ import annotations

import argparse
import json
import resource
import time
from pathlib import Path

import numpy as np
import pandas as pd

from fastde.abundance import run_abundance


def peak_rss_mb() -> float:
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0)


def make_fixture(n_samples: int, n_celltypes: int, cells_per_sample: int, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    samples = [f"S{i:03d}" for i in range(n_samples)]
    celltypes = ["RiskCells", "ResponderCells", "SubtypeACells"] + [f"Other{i}" for i in range(max(0, n_celltypes - 3))]
    labels = np.array(["response"] * (n_samples // 2) + ["non_response"] * (n_samples - n_samples // 2))
    rng.shuffle(labels)
    subtypes = np.resize(np.array(["control", "A", "B"]), n_samples)
    rng.shuffle(subtypes)
    counts = []
    risk_signal = []
    for i, sample in enumerate(samples):
        alpha = np.ones(len(celltypes))
        risk_strength = 0.1 + 0.75 * (i / max(n_samples - 1, 1))
        risk_signal.append(risk_strength)
        alpha[0] += 8 * risk_strength
        if labels[i] == "response":
            alpha[1] += 8
        if subtypes[i] == "A":
            alpha[2] += 8
        elif subtypes[i] == "B" and len(celltypes) > 3:
            alpha[3] += 8
        props = rng.dirichlet(alpha)
        counts.append(rng.multinomial(cells_per_sample, props))
    counts_df = pd.DataFrame(counts, index=samples, columns=celltypes)
    risk_signal = np.asarray(risk_signal)
    metadata = pd.DataFrame(
        {
            "sample_id": samples,
            "responder": labels,
            "subtype": subtypes,
            "OS_time": np.maximum(1.0, 120.0 - 90.0 * risk_signal + rng.normal(0, 4, n_samples)),
            "OS_event": (rng.random(n_samples) < 0.82).astype(int),
            "age": rng.normal(60, 8, n_samples),
            "sex": rng.choice(["F", "M"], n_samples),
        }
    )
    if metadata["OS_event"].sum() == 0:
        metadata.loc[0, "OS_event"] = 1
    return counts_df, metadata


def run_mode(
    mode: str,
    counts_path: Path,
    metadata_path: Path,
    output_dir: Path,
    *,
    epochs: int,
    learning_rate: float,
    max_instances_per_sample: int,
) -> dict[str, object]:
    start = time.perf_counter()
    kwargs = {
        "mode": mode,
        "counts": counts_path,
        "metadata": metadata_path,
        "covariates": "age,sex",
        "max_epochs": epochs,
        "learning_rate": learning_rate,
        "max_instances_per_sample": max_instances_per_sample,
        "output_dir": output_dir / mode,
    }
    if mode == "binary":
        kwargs.update({"label_col": "responder", "positive_label": "response", "negative_label": "non_response"})
    elif mode == "multiclass":
        kwargs.update({"label_col": "subtype", "reference_level": "control"})
    elif mode == "survival":
        kwargs.update({"time_col": "OS_time", "event_col": "OS_event", "survival_loss": "cox"})
    else:
        raise ValueError(mode)
    result = run_abundance(**kwargs)
    elapsed = time.perf_counter() - start
    top = result.results.sort_values(["padj", "pvalue"], na_position="last").iloc[0]["cell_type"]
    return {
        "mode": mode,
        "wall_seconds": elapsed,
        "peak_rss_mb": peak_rss_mb(),
        "top_cell_type": str(top),
        "metrics": result.metrics,
        "n_samples": int(result.manifest["n_samples"]),
        "n_celltypes": int(result.manifest["n_celltypes"]),
        "max_instances_per_sample": max_instances_per_sample,
        "epochs": epochs,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="tmp/fastde_abundance_mil_stress")
    parser.add_argument("--n-samples", type=int, default=72)
    parser.add_argument("--n-celltypes", type=int, default=8)
    parser.add_argument("--cells-per-sample", type=int, default=1000)
    parser.add_argument("--max-instances-per-sample", type=int, default=256)
    parser.add_argument("--epochs", type=int, default=40)
    parser.add_argument("--learning-rate", type=float, default=0.02)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    counts, metadata = make_fixture(args.n_samples, args.n_celltypes, args.cells_per_sample, args.seed)
    counts_path = out / "sample_by_celltype_counts.tsv"
    metadata_path = out / "sample_metadata.tsv"
    counts.to_csv(counts_path, sep="\t")
    metadata.to_csv(metadata_path, sep="\t", index=False)
    rows = []
    for mode in ["binary", "multiclass", "survival"]:
        rows.append(
            run_mode(
                mode,
                counts_path,
                metadata_path,
                out,
                epochs=args.epochs,
                learning_rate=args.learning_rate,
                max_instances_per_sample=args.max_instances_per_sample,
            )
        )
    flat_rows = []
    for row in rows:
        metrics = row["metrics"]
        flat = {key: value for key, value in row.items() if key != "metrics"}
        flat.update({f"metric_{key}": value for key, value in metrics.items() if isinstance(value, (int, float, str, bool))})
        flat_rows.append(flat)
    summary = pd.DataFrame(flat_rows)
    summary.to_csv(out / "stress_summary.tsv", sep="\t", index=False)
    (out / "stress_summary.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
