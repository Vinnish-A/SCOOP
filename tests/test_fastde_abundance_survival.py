from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from fastde.abundance import run_abundance
from fastde.cli import main


def _survival_fixture() -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(2)
    samples = [f"S{i}" for i in range(18)]
    counts = []
    risk_prop = []
    for i in range(18):
        p = 0.1 + 0.04 * i
        risk_prop.append(p)
        counts.append(rng.multinomial(400, [p, 1 - p]))
    counts_df = pd.DataFrame(counts, index=samples, columns=["RiskCells", "Other"])
    meta = pd.DataFrame({"sample_id": samples, "OS_time": 100 - 80 * np.array(risk_prop), "OS_event": [1] * 12 + [0] * 6})
    return counts_df, meta


def test_survival_mode_rejects_no_events(tmp_path) -> None:
    counts, meta = _survival_fixture()
    meta["OS_event"] = 0
    counts_path = tmp_path / "counts.tsv"
    meta_path = tmp_path / "metadata.tsv"
    counts.to_csv(counts_path, sep="\t")
    meta.to_csv(meta_path, sep="\t", index=False)
    with pytest.raises(ValueError, match="at least one observed event"):
        run_abundance(mode="survival", counts=counts_path, metadata=meta_path, time_col="OS_time", event_col="OS_event", output_dir=tmp_path / "out")


def test_cli_survival_writes_outputs(tmp_path) -> None:
    counts, meta = _survival_fixture()
    counts_path = tmp_path / "counts.tsv"
    meta_path = tmp_path / "metadata.tsv"
    counts.to_csv(counts_path, sep="\t")
    meta.to_csv(meta_path, sep="\t", index=False)
    rc = main(
        [
            "abundance",
            "--mode",
            "survival",
            "--counts",
            str(counts_path),
            "--metadata",
            str(meta_path),
            "--time-col",
            "OS_time",
            "--event-col",
            "OS_event",
            "--max-epochs",
            "200",
            "--learning-rate",
            "0.05",
            "--output-dir",
            str(tmp_path / "out"),
        ]
    )
    result = run_abundance(
        mode="survival",
        counts=counts_path,
        metadata=meta_path,
        time_col="OS_time",
        event_col="OS_event",
        max_epochs=200,
        learning_rate=0.05,
        output_dir=tmp_path / "out2",
    )
    assert rc == 0
    assert (tmp_path / "out" / "abundance_survival_results.tsv").exists()
    assert "concordance_index" in result.metrics
