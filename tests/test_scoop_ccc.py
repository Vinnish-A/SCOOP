from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml
from scipy import sparse

from scoop.commands import ccc
from scoop.results import result_payload

FAKE_FASTCCC = """
import argparse, anndata as ad
p = argparse.ArgumentParser(); p.add_argument("--input"); p.add_argument("--groupby"); p.add_argument("--output")
a = p.parse_args()
obs = ad.read_h5ad(a.input).obs
obs[a.groupby].astype(str).value_counts().rename_axis("group").reset_index().to_csv(a.output + "/groups.tsv", sep="\\t", index=False)
"""


def _run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[dict, Path, Path]:
    labels = ["T"] * 30 + ["B"] * 25 + ["Rare"] * 5 + ["Unresolved"] * 40
    status = ["accepted"] * 30 + ["labeled"] * 25 + ["accepted"] * 5 + ["unresolved"] * 40
    obs = pd.DataFrame({"cell_type_lvl3": labels, "cell_type_lvl2": labels, "annotation_status": status},
                       index=[f"c{i}" for i in range(100)])
    run_root = tmp_path / "runs/r1"
    (run_root / "artifacts").mkdir(parents=True)
    x = sparse.csr_matrix(np.ones((100, 3), dtype=np.float32))
    ad.AnnData(x, obs=obs, var=pd.DataFrame(index=["A", "B", "C"])).write_h5ad(run_root / "artifacts/adata_annotation_committed.h5ad")
    pd.DataFrame({"ligand": ["A"], "receptor": ["B|C"]}).to_csv(tmp_path / "lr.tsv", sep="\t", index=False)
    (tmp_path / "db.zip").write_bytes(b"")
    (tmp_path / "fake_fastccc.py").write_text(FAKE_FASTCCC, encoding="utf-8")
    config = run_root / "config/run.yaml"
    config.parent.mkdir(parents=True)
    config.write_text(yaml.safe_dump({"run": {"run_id": "r1"}, "ccc": {
        "enabled": True, "lr_resource": "lr.tsv", "fastccc_database": "db.zip", "min_cells_per_group": 20,
        "complex_validation": {"enabled": False},
        "fastccc_command_template": f"{sys.executable} {tmp_path / 'fake_fastccc.py'} --input {{input_h5ad}} "
                                    "--groupby {groupby} --output {output_dir}"}}), encoding="utf-8")
    monkeypatch.setattr(sys, "argv", ["ccc", "--config", str(config)])
    return ccc.main(), run_root, config


def test_ccc_tests_only_accepted_groups_with_enough_cells(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    result, run_root, config = _run(tmp_path, monkeypatch)

    groups = pd.read_csv(run_root / "ccc/fastccc/groups.tsv", sep="\t").set_index("group")["count"].to_dict()
    assert groups == {"T": 30, "B": 25}
    tested = ad.read_h5ad(run_root / "artifacts/adata_ccc_indexed.h5ad")
    assert tested.n_obs == 55 and list(tested.uns["ccc"]["dropped_groups"]) == ["Rare"]
    assert str((run_root / "ccc/fastccc/groups.tsv").resolve()) in result["outputs"]
    assert result["coverage"][0]["total"] == 100 and result["coverage"][0]["processed"] == 55
    assert result_payload(result, base=config.parent, route="managed")["execution"]["status"] == "completed"
