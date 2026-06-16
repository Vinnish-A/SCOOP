from __future__ import annotations

import numpy as np
import pandas as pd

from fastde.cli import main


def test_cli_binary_writes_outputs(tmp_path) -> None:
    samples = [f"S{i}" for i in range(12)]
    labels = ["yes"] * 6 + ["no"] * 6
    counts = pd.DataFrame(
        [[80, 20] if label == "yes" else [20, 80] for label in labels],
        index=samples,
        columns=["A", "B"],
    )
    meta = pd.DataFrame({"sample_id": samples, "responder": labels, "age": np.arange(12)})
    counts_path = tmp_path / "counts.tsv"
    meta_path = tmp_path / "metadata.tsv"
    counts.to_csv(counts_path, sep="\t")
    meta.to_csv(meta_path, sep="\t", index=False)
    out = tmp_path / "out"
    rc = main(
        [
            "abundance",
            "--mode",
            "binary",
            "--counts",
            str(counts_path),
            "--metadata",
            str(meta_path),
            "--label-col",
            "responder",
            "--positive-label",
            "yes",
            "--negative-label",
            "no",
            "--covariates",
            "age",
            "--max-epochs",
            "100",
            "--learning-rate",
            "0.05",
            "--output-dir",
            str(out),
        ]
    )
    assert rc == 0
    assert (out / "abundance_binary_results.tsv").exists()
    assert (out / "abundance_manifest.json").exists()
