from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from fastde.abundance import run_abundance


def _fixture(n_samples: int = 60, n_celltypes: int = 8, seed: int = 0) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    celltypes = [f"CT{i}" for i in range(n_celltypes)]
    samples = [f"S{i:03d}" for i in range(n_samples)]
    labels = np.array(["response"] * (n_samples // 2) + ["non_response"] * (n_samples - n_samples // 2))
    rng.shuffle(labels)
    subtype = np.array(["control", "A", "B"] * (n_samples // 3) + ["control"] * (n_samples % 3))[:n_samples]
    rng.shuffle(subtype)
    base = np.ones(n_celltypes)
    counts = []
    for i in range(n_samples):
        alpha = base.copy()
        if labels[i] == "response":
            alpha[0] += 5
        if subtype[i] == "A":
            alpha[1] += 5
        if subtype[i] == "B":
            alpha[2] += 5
        props = rng.dirichlet(alpha)
        counts.append(rng.multinomial(1000, props))
    counts_df = pd.DataFrame(counts, index=samples, columns=celltypes)
    risk = counts_df["CT3"] / counts_df.sum(axis=1)
    os_time = np.maximum(1, 100 - 60 * risk.to_numpy() + rng.normal(0, 5, n_samples))
    os_event = (rng.random(n_samples) < 0.75).astype(int)
    meta = pd.DataFrame(
        {
            "sample_id": samples,
            "responder": labels,
            "subtype": subtype,
            "OS_time": os_time,
            "OS_event": os_event,
            "age": rng.normal(60, 8, n_samples),
            "sex": rng.choice(["F", "M"], n_samples),
        }
    )
    return counts_df, meta


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="tmp/fastde_abundance_benchmark")
    parser.add_argument("--n-samples", type=int, default=60)
    parser.add_argument("--n-celltypes", type=int, default=8)
    args = parser.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    counts, meta = _fixture(args.n_samples, args.n_celltypes)
    counts_path = out / "sample_by_celltype_counts.tsv"
    meta_path = out / "sample_metadata.tsv"
    counts.to_csv(counts_path, sep="\t")
    meta.to_csv(meta_path, sep="\t", index=False)
    run_abundance(
        mode="binary",
        counts=counts_path,
        metadata=meta_path,
        label_col="responder",
        positive_label="response",
        negative_label="non_response",
        covariates="age,sex",
        output_dir=out / "binary",
    )
    run_abundance(
        mode="multiclass",
        counts=counts_path,
        metadata=meta_path,
        label_col="subtype",
        reference_level="control",
        covariates="age,sex",
        output_dir=out / "multiclass",
    )
    run_abundance(
        mode="survival",
        counts=counts_path,
        metadata=meta_path,
        time_col="OS_time",
        event_col="OS_event",
        covariates="age,sex",
        output_dir=out / "survival",
    )
    print(f"wrote benchmark outputs to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
