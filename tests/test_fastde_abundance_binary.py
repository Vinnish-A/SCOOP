from __future__ import annotations

import numpy as np
import pandas as pd

from fastde.abundance import run_abundance


def _binary_fixture() -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(1)
    samples = [f"S{i}" for i in range(20)]
    labels = np.array(["response"] * 10 + ["non_response"] * 10)
    counts = []
    for label in labels:
        alpha = np.array([8, 1, 1]) if label == "response" else np.array([1, 8, 1])
        counts.append(rng.multinomial(500, rng.dirichlet(alpha)))
    return pd.DataFrame(counts, index=samples, columns=["T_enriched", "B_enriched", "Other"]), pd.DataFrame({"sample_id": samples, "responder": labels})


def test_binary_mil_loss_decreases(tmp_path) -> None:
    counts, meta = _binary_fixture()
    counts_path = tmp_path / "counts.tsv"
    meta_path = tmp_path / "metadata.tsv"
    counts.to_csv(counts_path, sep="\t")
    meta.to_csv(meta_path, sep="\t", index=False)
    result = run_abundance(
        mode="binary",
        counts=counts_path,
        metadata=meta_path,
        label_col="responder",
        positive_label="response",
        negative_label="non_response",
        max_epochs=80,
        learning_rate=0.05,
        output_dir=tmp_path / "loss",
    )
    history = result.manifest["model"]["history"]
    assert history[-1]["loss"] < history[0]["loss"]


def test_binary_mode_recovers_enriched_cell_type(tmp_path) -> None:
    counts, meta = _binary_fixture()
    counts_path = tmp_path / "counts.tsv"
    meta_path = tmp_path / "metadata.tsv"
    counts.to_csv(counts_path, sep="\t")
    meta.to_csv(meta_path, sep="\t", index=False)
    result = run_abundance(
        mode="binary",
        counts=counts_path,
        metadata=meta_path,
        label_col="responder",
        positive_label="response",
        negative_label="non_response",
        max_epochs=300,
        learning_rate=0.05,
        output_dir=tmp_path / "out",
    )
    top = result.results.sort_values("pvalue").iloc[0]["cell_type"]
    assert top in {"T_enriched", "B_enriched"}
    assert result.metrics["auc"] >= 0.9
