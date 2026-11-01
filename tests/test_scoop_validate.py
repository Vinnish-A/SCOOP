from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from scipy import sparse

from scoop.commands import validate
from scoop.results import result_payload


def _project(tmp_path: Path, *, obs_names=None, var_names=None, counts=None) -> Path:
    x = np.arange(30, dtype=np.float32).reshape(6, 5) % 4 if counts is None else counts
    obs = pd.DataFrame({"sample_id": ["s1"] * 6}, index=obs_names or [f"c{i}" for i in range(6)])
    var = pd.DataFrame(index=var_names or [f"G{i}" for i in range(5)])
    ad.AnnData(sparse.csr_matrix(x), obs=obs, var=var, layers={"counts": sparse.csr_matrix(x)}).write_h5ad(tmp_path / "in.h5ad")
    config = tmp_path / "runs/r1/config/run.yaml"
    config.parent.mkdir(parents=True)
    config.write_text(yaml.safe_dump({"run": {"run_id": "r1", "input_h5ad": "in.h5ad"}}), encoding="utf-8")
    return config


def _run(monkeypatch: pytest.MonkeyPatch, config: Path) -> dict:
    monkeypatch.setattr(sys, "argv", ["validate", "--config", str(config)])
    return validate.main()


def test_validate_registers_output_and_returns_declared_artifacts(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config = _project(tmp_path)

    result = _run(monkeypatch, config)

    output = tmp_path / "runs/r1/artifacts/validated_adata.h5ad"
    record = ad.read_h5ad(output).uns["file_registry"]["artifacts"]["validated_adata"]
    assert Path(record["path"]) == output
    assert result["outputs"] == [str(output.resolve())]
    assert result["inputs"] == [str((tmp_path / "in.h5ad").resolve())]
    assessed = result_payload(result, base=config.parent, route="managed")
    assert assessed["validation"]["status"] == "passed"


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"obs_names": ["c1", "c1", "c2", "c3", "c4", "c5"]}, "obs_names are not unique"),
        ({"var_names": ["MT-CO1", "GAPDH", "GAPDH", "ACTB", "B2M"]}, "var_names are not unique"),
    ],
)
def test_validate_rejects_duplicate_names(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, kwargs: dict, message: str) -> None:
    config = _project(tmp_path, **kwargs)
    with pytest.raises(ValueError, match=message):
        _run(monkeypatch, config)
    assert not (tmp_path / "runs/r1/artifacts/validated_adata.h5ad").exists()


@pytest.mark.parametrize("bad", [2.5, -1.0, np.nan, np.inf])
def test_validate_rejects_non_integer_counts(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, bad: float) -> None:
    counts = np.ones((6, 5), dtype=np.float32)
    counts[0, 0] = bad
    config = _project(tmp_path, counts=counts)
    with pytest.raises(ValueError, match="non-negative integer counts"):
        _run(monkeypatch, config)
