from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import scanpy as sc
from anndata import AnnData

from scsp_agent_sop.core import leiden_sweep
from scsp_agent_sop.core_runner import run_core_pipeline
from scsp_agent_sop.storage import init_file_registry


def test_core_runner_writes_manifest_with_planned_fallback(tmp_path, monkeypatch):
    adata = AnnData(np.ones((5, 4)))
    init_file_registry(adata, "test")

    def fake_legacy(adata, cfg, run_root):
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
    manifest = json.loads((tmp_path / "02_core" / "fastcore" / "fastcore_manifest.json").read_text())
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
