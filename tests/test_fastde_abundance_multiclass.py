from __future__ import annotations

import numpy as np
import pandas as pd

from fastde.abundance import run_abundance


def test_multiclass_mode_runs(tmp_path) -> None:
    rng = np.random.default_rng(3)
    samples = [f"S{i}" for i in range(24)]
    classes = np.array(["control", "A", "B"] * 8)
    counts = []
    for cls in classes:
        alpha = np.array([6, 1, 1]) if cls == "control" else np.array([1, 6, 1]) if cls == "A" else np.array([1, 1, 6])
        counts.append(rng.multinomial(500, rng.dirichlet(alpha)))
    counts_df = pd.DataFrame(counts, index=samples, columns=["C0", "C1", "C2"])
    meta = pd.DataFrame({"sample_id": samples, "subtype": classes})
    counts_path = tmp_path / "counts.tsv"
    meta_path = tmp_path / "metadata.tsv"
    counts_df.to_csv(counts_path, sep="\t")
    meta.to_csv(meta_path, sep="\t", index=False)
    result = run_abundance(
        mode="multiclass",
        counts=counts_path,
        metadata=meta_path,
        label_col="subtype",
        reference_level="control",
        max_epochs=300,
        learning_rate=0.05,
        output_dir=tmp_path / "out",
    )
    assert result.results["class_or_contrast"].nunique() == 2
    assert result.metrics["accuracy"] >= 0.8
