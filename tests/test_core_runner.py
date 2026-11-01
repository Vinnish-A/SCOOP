from __future__ import annotations

import json

import numpy as np
import scanpy as sc
from anndata import AnnData

from scsp_agent_sop.core import leiden_sweep
from scsp_agent_sop.core_runner import CORE_OUTPUTS, run_core_pipeline
from scsp_agent_sop.storage import init_file_registry


def _deliver_core_outputs(adata):
    """What a real backend writes; the runner refuses to report success without it."""
    adata.layers["counts"] = adata.X.copy()
    adata.layers["log1p_norm"] = adata.X.copy()
    for key in CORE_OUTPUTS["obsm"]:
        adata.obsm[key] = np.ones((adata.n_obs, 2))
    adata.obsp["connectivities_identity"] = np.eye(adata.n_obs)
    adata.obs["cluster_identity"] = "0"
    adata.var["highly_variable_biology"] = True


def test_core_runner_writes_manifest_with_planned_fallback(tmp_path, monkeypatch):
    adata = AnnData(np.ones((5, 4)))
    init_file_registry(adata, "test")

    def fake_legacy(adata, cfg, run_root):
        _deliver_core_outputs(adata)
        adata.obs["cluster_identity"] = ["0", "0", "1", "1", "1"]
        return {
            "backend": "scanpy_legacy",
            "batch_correction_method": "harmony2",
            "harmony2_used": False,
            "batch_keys": [],
            "n_obs": adata.n_obs,
            "n_vars": adata.n_vars,
            "n_clusters": 2,
            "timings": {"fake": 0.1},
            "quality": {"accepted": True},
            "artifacts": {},
        }

    monkeypatch.setattr("scsp_agent_sop.core_runner.run_scanpy_legacy_core", fake_legacy)
    cfg = {
        "run": {"run_id": "test"},
        "core": {
            "engine": "fastcore",
            "fallback_engine": "scanpy_legacy",
            "fastcore": {"fallback_backend": "scanpy_legacy", "enable_fastcore_cpu_backend": False},
        },
    }
    result = run_core_pipeline(adata, cfg, tmp_path)
    assert result["backend"] == "scanpy_legacy"
    assert result["fallback_used"] is True
    manifest = json.loads((tmp_path / "core" / "fastcore" / "fastcore_manifest.json").read_text())
    assert manifest["backend"] == "scanpy_legacy"
    assert manifest["fallback_used"] is True
    assert "fastcore_manifest" in adata.uns["file_registry"]["artifacts"]
    assert (tmp_path / "logs" / "decision_log.jsonl").exists()


def test_leiden_sweep_uses_coarse_to_fine_search() -> None:
    rng = np.random.default_rng(0)
    adata = AnnData(rng.normal(size=(36, 6)))
    sc.pp.neighbors(adata, n_neighbors=6, random_state=0)
    adata.obsp["connectivities_identity"] = adata.obsp["connectivities"].copy()

    resolutions = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5]
    seeds = [0, 1, 2, 3, 4]
    sweep, stability = leiden_sweep(
        adata,
        resolutions=resolutions,
        seeds=seeds,
        search_config={
            "strategy": "coarse_to_fine",
            "coarse_resolutions": [0.25, 0.75, 1.25, 1.5],
            "coarse_seeds": [0, 1],
            "refine_window": 0.25,
            "refine_step": 0.125,
            "min_seed_ari": 0.0,
        },
    )

    assert "cluster_identity" in adata.obs
    assert len(sweep) < len(resolutions) * len(seeds)
    assert {"coarse", "refine"}.issubset(set(sweep["phase"]))
    assert stability["chosen"].sum() == 1
    assert set(stability["search_strategy"]) == {"coarse_to_fine"}


def test_core_runner_filters_qc_failures_before_backend(tmp_path, monkeypatch):
    adata = AnnData(np.ones((5, 4)))
    adata.obs["final_use"] = [True, False, True, True, False]
    init_file_registry(adata, "test")

    def backend(adata, cfg, run_root):
        assert adata.n_obs == 3
        assert adata.obs.final_use.all()
        _deliver_core_outputs(adata)
        return {"backend": "scanpy_legacy", "n_obs": 3, "timings": {}, "artifacts": {}}

    monkeypatch.setattr("scsp_agent_sop.core_runner.run_scanpy_legacy_core", backend)
    run_core_pipeline(adata, {"core": {"engine": "scanpy_legacy"}}, tmp_path)
    manifest = json.loads((tmp_path / "core/fastcore/fastcore_manifest.json").read_text())
    assert manifest["qc_input_observations"] == 5
    assert manifest["n_obs"] == 3


def test_core_runner_refuses_success_without_its_output_contract(tmp_path, monkeypatch):
    import pytest

    adata = AnnData(np.ones((5, 4)))
    init_file_registry(adata, "test")

    def incomplete(adata, cfg, run_root):
        _deliver_core_outputs(adata)
        del adata.layers["log1p_norm"]
        return {"backend": "scanpy_legacy", "timings": {}, "artifacts": {}}

    monkeypatch.setattr("scsp_agent_sop.core_runner.run_scanpy_legacy_core", incomplete)
    with pytest.raises(RuntimeError, match=r"missing layers\['log1p_norm'\]"):
        run_core_pipeline(adata, {"core": {"engine": "scanpy_legacy"}}, tmp_path)
    assert not (tmp_path / "core/fastcore/fastcore_manifest.json").exists()


def _sweep_object():
    rng = np.random.default_rng(0)
    centres = rng.normal(scale=6, size=(4, 6))
    adata = AnnData(np.vstack([centre + rng.normal(size=(15, 6)) for centre in centres]))
    sc.pp.neighbors(adata, n_neighbors=8, random_state=0)
    adata.obsp["connectivities_identity"] = adata.obsp["connectivities"].copy()
    return adata


def test_leiden_selection_policy_is_recorded_and_sweep_labels_are_kept():
    resolutions = [0.1, 0.5, 1.0, 2.0]
    chosen = {}
    for policy in ("coarsest_stable", "finest_stable"):
        adata = _sweep_object()
        _, stability = leiden_sweep(adata, resolutions=resolutions, seeds=[0, 1],
                                    search_config={"strategy": "full", "min_seed_ari": 0.0, "selection": policy})
        assert stability.eligible.all() and set(stability.selection_policy) == {policy}
        record = adata.uns["leiden_selection"]["cluster_identity"]
        chosen[policy] = record["chosen_resolution"]
        assert record["policy"] == policy and record["rule_satisfied"]
        assert [row["resolution"] for row in json.loads(record["candidates_json"])] == resolutions
        labels = adata.obsm["cluster_identity_sweep_labels"]
        assert list(labels.columns) == [f"r{res}" for res in resolutions]
        assert (labels[f"r{chosen[policy]}"].astype(str) == adata.obs.cluster_identity.astype(str)).all()
        assert not any(column.startswith("_tmp_leiden_") for column in adata.obs)
    assert chosen == {"coarsest_stable": 0.1, "finest_stable": 2.0}


def test_leiden_finest_policy_refines_around_the_finest_stable_coarse_resolution():
    adata = _sweep_object()
    sweep, _ = leiden_sweep(adata, resolutions=[0.25, 0.5, 0.75, 1.0, 1.25, 1.5], seeds=[0, 1, 2],
                            search_config={"strategy": "coarse_to_fine", "coarse_resolutions": [0.25, 1.5],
                                           "coarse_seeds": [0, 1], "refine_window": 0.25, "refine_step": 0.25,
                                           "min_seed_ari": 0.0, "selection": "finest_stable"})
    assert set(sweep.loc[sweep.phase == "refine", "resolution"]) == {1.25, 1.5}
    assert adata.uns["leiden_selection"]["cluster_identity"]["chosen_resolution"] == 1.5


def test_core_manifest_reports_a_harmony_import_fallback_for_review(tmp_path, monkeypatch):
    from scipy import sparse

    def unavailable(*args, **kwargs):
        raise ImportError("Harmony 2.0 is not available.")

    monkeypatch.setattr("scsp_agent_sop.core.run_harmony2", unavailable)
    rng = np.random.default_rng(0)
    counts = sparse.csr_matrix(rng.poisson(1.2, size=(40, 60)).astype(np.float32))
    adata = AnnData(counts.copy())
    adata.var_names = [f"Gene{i}" for i in range(adata.n_vars)]
    adata.obs_names = [f"Cell{i}" for i in range(adata.n_obs)]
    adata.obs["dataset"] = ["d1", "d2"] * 20
    adata.layers["counts"] = counts.copy()
    init_file_registry(adata, "test")
    cfg = {"run": {"obs_keys": {"batch": "dataset"}},
           "core": {"engine": "fastcore", "fastcore": {"allowed_backends": ["fastcore_cpu"],
                                                       "omicverse": {"n_hvgs": 20, "n_pcs": 5,
                                                                     "leiden": {"resolutions": [0.6]}}}}}
    result = run_core_pipeline(adata, cfg, tmp_path)
    assert result["harmony2_used"] is False and result["batch_keys"] == []
    manifest = json.loads((tmp_path / "core/fastcore/fastcore_manifest.json").read_text())
    record = manifest["batch_correction"]
    assert record["status"] == "fallback_uncorrected" and record["requested_keys"] == ["dataset"]
    assert record["key_source"] == "run.obs_keys.batch" and "not available" in record["reason"]
    np.testing.assert_array_equal(adata.obsm["X_pca_harmony_identity"], adata.obsm["X_pca_identity_prebatch"])
    decision = json.loads((tmp_path / "logs/decision_log.jsonl").read_text().splitlines()[-1])
    assert decision["human_review_required"] and "uncorrected PCA" in decision["review_reason"]
    assert manifest["clustering"]["policy"] == "coarsest_stable" and manifest["clustering"]["candidates"]
