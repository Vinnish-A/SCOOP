import numpy as np
import pytest
from anndata import AnnData, read_h5ad

from scsp_agent_sop.programs import run_fast_consensus_nmf, run_nmf_once, consensus_stability, NMFRunResult, choose_k
import pandas as pd


def test_replicate_seeds_have_independent_starts_and_report_budget():
    x = np.random.default_rng(8).gamma(2.0, 1.0, (40, 30))
    with pytest.warns(UserWarning):
        first = run_nmf_once(x, 3, 0, max_iter=1)
    with pytest.warns(UserWarning):
        second = run_nmf_once(x, 3, 1, max_iter=1)
    assert not np.allclose(first.H, second.H)
    assert first.n_iter == 1 and not first.converged


def test_native_programs_deliver_declared_consensus_and_convergence(tmp_path):
    rng = np.random.default_rng(12)
    a = AnnData(rng.gamma(2.0, 1.0, (50, 24)))
    a.layers["log1p_norm"] = a.X.copy()
    summary, weights, usage = run_fast_consensus_nmf(a, k_grid=[2, 3], seeds=[0, 1, 2], max_iter=30, workers=2)
    assert summary.method.eq("sklearn_multistart_nmf").all()
    assert summary.init.eq("random").all()
    assert summary.converged_fraction.between(0, 1).all()
    assert summary.median_iterations.le(30).all()
    assert summary.selected.sum() == 1
    k = int(summary.loc[summary.selected, "k"].iloc[0])
    assert weights.shape == (k, 25)
    assert usage.shape == (50, k + 1)
    assert a.uns["nmf_method"]["output"] == "consensus_components_and_fixed_H_usage_refit"
    a.write_h5ad(tmp_path / "programs.h5ad")
    assert read_h5ad(tmp_path / "programs.h5ad").obsm["X_nmf_usage"].shape == (50, k)
    assert summary.min_component_stability.le(summary.median_stability).all()
    assert a.uns["nmf_method"]["discovery_features"] == 24


def _run(H, seed):
    H = np.asarray(H, dtype=float)
    return NMFRunResult(k=len(H), seed=seed, error=1., W=np.ones((4, len(H))), H=H, n_iter=2, converged=True)


def test_same_seed_duplicates_and_singletons_cannot_establish_stability():
    first = _run([[1, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0]], 0)
    second = _run([[1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], 1)
    stability, consensus = consensus_stability([first, second], 3)
    assert stability == 0 and consensus.shape == (3, 4)
    assert np.isnan(consensus_stability([first], 3)[0])
    with pytest.raises(ValueError, match="independently seeded"):
        consensus_stability([first, first], 3)


def test_component_permutations_match_across_independent_seeds():
    a = _run(np.eye(3), 0)
    b = _run(np.eye(3)[[2, 0, 1]], 1)
    score, components = consensus_stability([a, b], 3)
    assert score == 1
    np.testing.assert_allclose(components, np.eye(3))


def test_k_selection_does_not_hide_a_weak_component_under_high_median():
    summary = pd.DataFrame({"k": [2, 3], "median_stability": [.99, .95],
                            "min_component_stability": [.1, .9], "median_error": [1., .8]})
    assert choose_k(summary, .7) == 3


def test_programs_record_an_unhonoured_measurement_and_refuse_a_guessed_expression_batch(tmp_path, monkeypatch):
    import sys

    import yaml
    from scoop.commands import programs
    from scsp_agent_sop.measurement import FORMAT, FeatureMeasurement
    from scsp_agent_sop.storage import init_file_registry, register_file

    rng = np.random.default_rng(12)
    a = AnnData(rng.gamma(2.0, 1.0, (40, 20)), obs=pd.DataFrame({"sample_id": ["s1", "s2"] * 20},
                                                               index=[f"c{i}" for i in range(40)]))
    a.layers["log1p_norm"] = a.X.copy()
    a.layers["counts"] = np.rint(a.X * 3)
    init_file_registry(a, "r1")
    fm = FeatureMeasurement.from_var("sample_id", ["s1", "s2"], a.var_names, np.ones((2, 20), bool))
    register_file(a, key="feature_measurement", path=fm.save(tmp_path / "fm.npz"), category="artifacts", schema=FORMAT)
    run = tmp_path / "runs" / "r1"
    (run / "config").mkdir(parents=True)
    (run / "artifacts").mkdir()
    a.write_h5ad(run / "artifacts" / "adata_core.h5ad")
    config = run / "config" / "run.yaml"
    nmf = {"k_grid": [2], "seeds": [0, 1], "max_iter": 30, "workers": 1}
    config.write_text(yaml.safe_dump({"programs": nmf}))
    monkeypatch.setattr(sys, "argv", ["scoop programs", "--config", str(config)])
    receipt = programs.main()
    assert receipt["measurement"] == {"status": "not_honoured", "path": str(tmp_path / "fm.npz")}
    assert any("feature measurement" in item for item in receipt["limitations"])

    config.write_text(yaml.safe_dump({"programs": {**nmf, "batch_correction": "harmony_moe"},
                                      "keys": {"batch_candidates": ["sample_id"]}}))
    with pytest.raises(ValueError, match="declared run.obs_keys.batch"):
        programs.main()
