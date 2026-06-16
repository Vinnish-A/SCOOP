from __future__ import annotations

import numpy as np
import pandas as pd

from fastde.abundance import run_abundance
from fastde.abundance_data import AbundanceTable
from fastde.abundance_design import build_feature_design
from fastde.abundance_train import AbundanceTrainer, AbundanceTrainingConfig
from fastde.abundance_loss import binary_bce_with_logits


def _binary_fixture() -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(1)
    samples = [f"S{i}" for i in range(20)]
    labels = np.array(["response"] * 10 + ["non_response"] * 10)
    counts = []
    for label in labels:
        alpha = np.array([8, 1, 1]) if label == "response" else np.array([1, 8, 1])
        counts.append(rng.multinomial(500, rng.dirichlet(alpha)))
    return pd.DataFrame(counts, index=samples, columns=["T_enriched", "B_enriched", "Other"]), pd.DataFrame({"sample_id": samples, "responder": labels})


def test_binary_loss_decreases() -> None:
    counts, meta = _binary_fixture()
    table = AbundanceTable(counts, counts.div(counts.sum(axis=1), axis=0), meta.set_index("sample_id"))
    design = build_feature_design(table)
    y = (meta["responder"] == "response").astype(float).to_numpy()
    model = AbundanceTrainer("binary", AbundanceTrainingConfig(max_epochs=200, learning_rate=0.05)).fit(design.features, y=y)
    before = binary_bce_with_logits(np.zeros_like(y), y)
    after = binary_bce_with_logits(model.predict_score(design.features.to_numpy()), y)
    assert after < before


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
