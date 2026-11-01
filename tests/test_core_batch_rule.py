"""Core Harmony batch rule (declared run.obs_keys.batch, else candidates) and what Harmony reports."""
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scsp_agent_sop.core import core_batch_keys, correct_identity_pca, run_harmony2


def _obs():
    return pd.DataFrame({"sample_id": ["s1", "s2", "s3", "s4"] * 15, "dataset": ["d1"] * 30 + ["d2"] * 30,
                         "chemistry": ["v3"] * 60}, index=[f"c{i}" for i in range(60)])


def test_a_declared_batch_is_used_exactly_and_candidates_are_recorded_as_such():
    candidates = {"keys": {"batch_candidates": ["sample_id", "library_id", "chemistry"]}}
    assert core_batch_keys(candidates, _obs()) == (["sample_id", "chemistry"], "keys.batch_candidates")
    declared = {**candidates, "run": {"obs_keys": {"batch": "dataset"}}}
    assert core_batch_keys(declared, _obs()) == (["dataset"], "run.obs_keys.batch")
    assert core_batch_keys({}, _obs()) == (["sample_id"], "sample_key_default")
    with pytest.raises(KeyError, match="run.obs_keys"):
        core_batch_keys({"run": {"obs_keys": {"batch": "lane"}}}, _obs())


def _pca_object():
    rng = np.random.default_rng(0)
    adata = AnnData(np.zeros((60, 2)), obs=_obs())
    adata.obsm["X_pca_identity_prebatch"] = rng.normal(size=(60, 5)) + (adata.obs.dataset == "d2").to_numpy()[:, None]
    return adata


def test_harmony_reports_the_keys_it_used_its_iterations_and_objective():
    adata = _pca_object()
    record = run_harmony2(adata, batch_keys=["dataset", "chemistry"], max_iter_harmony=5)
    assert record["status"] == "corrected" and record["keys"] == ["dataset"]
    assert record["dropped_keys"] == {"chemistry": "single_level"}
    assert 1 <= record["n_iter"] <= 5 and len(record["objective_harmony"]) == record["n_iter"] + 1
    assert record["package"].startswith("harmonypy 2.")
    assert not np.allclose(adata.obsm["X_pca_harmony_identity"], adata.obsm["X_pca_identity_prebatch"])


def test_harmony_with_no_two_level_key_says_it_did_not_run():
    adata = _pca_object()
    record = correct_identity_pca(adata, method="harmony2", batch_keys=["chemistry"], key_source="keys.batch_candidates")
    assert record["status"] == "not_run" and record["keys"] == [] and "two levels" in record["reason"]
    assert record["requested_keys"] == ["chemistry"] and record["key_source"] == "keys.batch_candidates"
    np.testing.assert_array_equal(adata.obsm["X_pca_harmony_identity"], adata.obsm["X_pca_identity_prebatch"])
    assert correct_identity_pca(adata, method="none", batch_keys=["dataset"], key_source="x")["status"] == "not_requested"
