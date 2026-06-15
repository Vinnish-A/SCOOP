from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from fastcnmf.quality import evaluate_benchmark_gate


def _write_outputs(base: Path, run_name: str, k: int, order: list[str]) -> None:
    base.mkdir(parents=True, exist_ok=True)
    spectra = pd.DataFrame(
        {
            "gene_a": [1.0, 0.0],
            "gene_b": [0.0, 1.0],
            "gene_c": [0.5, 0.2],
        },
        index=["1", "2"],
    ).loc[order]
    usage = pd.DataFrame(
        {
            "1": [0.9, 0.8, 0.1, 0.2],
            "2": [0.1, 0.2, 0.9, 0.8],
        },
        index=[f"cell{i}" for i in range(4)],
    )
    usage = usage.loc[:, order]
    spectra.to_csv(base / f"{run_name}.spectra.k_{k}.dt_0_5.consensus.txt", sep="\t")
    usage.to_csv(base / f"{run_name}.usages.k_{k}.dt_0_5.consensus.txt", sep="\t")


def test_benchmark_gate_matches_permuted_programs(tmp_path: Path) -> None:
    summary = {
        "comparison": {
            "serial_factorize_seconds": 100.0,
            "parallel_factorize_seconds": 40.0,
        }
    }
    summary_path = tmp_path / "summary.json"
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    ref = tmp_path / "ref"
    cand = tmp_path / "cand"
    _write_outputs(ref, "run", 2, ["1", "2"])
    _write_outputs(cand, "run", 2, ["2", "1"])

    result = evaluate_benchmark_gate(
        summary_json=summary_path,
        reference_dir=ref,
        candidate_dir=cand,
        run_name="run",
        k_values=(2,),
    )
    assert result.passed
    assert result.time_saved_fraction == 0.6
    assert result.consistency[0].overall_consistency > 0.99

