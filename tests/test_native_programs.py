import numpy as np
import pytest
from anndata import AnnData, read_h5ad

from scsp_agent_sop.programs import run_fast_consensus_nmf, run_nmf_once


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
