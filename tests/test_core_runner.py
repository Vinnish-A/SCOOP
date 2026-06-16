from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from anndata import AnnData

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
