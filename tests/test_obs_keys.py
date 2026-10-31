from __future__ import annotations

import sys

import numpy as np
import pandas as pd
import pytest
import yaml
from anndata import AnnData

from scsp_agent_sop.config import resolve_obs_keys


def _obs(**columns):
    return pd.DataFrame(columns or {"sample_id": ["s1"]})


def test_obs_keys_fall_back_to_legacy_keys_and_record_sources():
    cfg = {"keys": {"sample": "lib", "batch_candidates": ["library_id", "chemistry"]}}
    keys = resolve_obs_keys(cfg, _obs(lib=["a"], chemistry=["v3"]))
    assert (keys.sample, keys.donor, keys.batch, keys.study) == ("lib", "donor_id", "chemistry", "study_id")
    assert keys.sources == {"sample": "keys.sample", "donor": "default", "batch": "keys.batch_candidates", "study": "default"}
    cfg["run"] = {"obs_keys": {"sample": "sample_id", "batch": "library_id", "donor": None}}
    keys = resolve_obs_keys(cfg, _obs(sample_id=["a"], library_id=["L1"]),
                            overrides={"sample": "annotation.decision.sample_key"})
    assert (keys.sample, keys.batch) == ("sample_id", "library_id")
    assert keys.record()["batch"] == {"key": "library_id", "source": "run.obs_keys.batch"}
    cfg["annotation"] = {"decision": {"sample_key": "lib"}}
    keys = resolve_obs_keys(cfg, _obs(lib=["a"], library_id=["L1"]), overrides={"sample": "annotation.decision.sample_key"})
    assert keys.record()["sample"] == {"key": "lib", "source": "annotation.decision.sample_key"}


def test_declared_or_required_obs_keys_must_exist():
    with pytest.raises(KeyError, match=r"batch='lane' \(run.obs_keys.batch\)"):
        resolve_obs_keys({"run": {"obs_keys": {"batch": "lane"}}}, _obs(sample_id=["a"]))
    with pytest.raises(KeyError, match=r"batch=None \(unset\)"):
        resolve_obs_keys({}, _obs(sample_id=["a"]), required=("batch",))
    # legacy and default names may be absent, as before R10
    assert resolve_obs_keys({}, _obs(other=["a"])).donor == "donor_id"


def _run(tmp_path, monkeypatch, cfg):
    run = tmp_path / "runs" / "test"
    (run / "config").mkdir(parents=True)
    (run / "artifacts").mkdir()
    config = run / "config" / "run.yaml"
    config.write_text(yaml.safe_dump(cfg))
    monkeypatch.chdir(tmp_path)

    def call(module, *args):
        monkeypatch.setattr(sys, "argv", ["scoop", "--config", str(config), *args])
        return module.main()
    return run, call


def test_de_prepare_aggregates_accepted_cells_with_design_columns(tmp_path, monkeypatch):
    from scoop.commands import de_prepare

    n = 8
    obs = pd.DataFrame({"lib": ["s1", "s1", "s1", "s1", "s2", "s2", "s2", "s2"], "condition": ["a"] * 4 + ["b"] * 4,
                        "cell_type_lvl3": ["T"] * n, "lane": ["L1"] * 4 + ["L2"] * 4, "chemistry": ["v3"] * n,
                        "annotation_status": ["accepted", "accepted", "accepted", "review_required"] + ["accepted"] * 3 + ["unassigned"]},
                       index=[f"c{i}" for i in range(n)])
    adata = AnnData(np.ones((n, 2)), obs=obs, var=pd.DataFrame(index=["G1", "G2"]))
    adata.layers["counts"] = adata.X.copy()
    run, call = _run(tmp_path, monkeypatch, {"run": {"obs_keys": {"sample": "lib", "batch": "lane"}},
                                             "keys": {"batch_candidates": ["sample_id", "chemistry"]},
                                             "de": {"min_cells_per_sample_celltype": 1}})
    adata.write_h5ad(run / "artifacts" / "final_adata.h5ad")
    call(de_prepare)
    meta = pd.read_csv(run / "differential-expression" / "pseudobulk" / "T" / "metadata.tsv", sep="\t")
    assert meta.n_cells.tolist() == [3, 3]
    assert meta.lane.tolist() == ["L1", "L2"] and meta.chemistry.tolist() == ["v3", "v3"]
    counts = pd.read_csv(run / "differential-expression" / "pseudobulk" / "T" / "counts.tsv", sep="\t")
    assert counts.G1.tolist() == [3.0, 3.0]

    run.joinpath("config", "run.yaml").write_text(yaml.safe_dump({"run": {"obs_keys": {"sample": "library"}}}))
    with pytest.raises(KeyError, match="sample='library'"):
        call(de_prepare)
