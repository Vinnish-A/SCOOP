from __future__ import annotations

import numpy as np
from anndata import AnnData
from scipy import sparse

from scsp_agent_sop.core_runner import run_core_pipeline
from scsp_agent_sop.storage import init_file_registry


def test_vendored_omicverse_cpu_backend_writes_stable_core_fields(tmp_path):
    rng = np.random.default_rng(0)
    counts = sparse.csr_matrix(rng.poisson(1.2, size=(40, 80)).astype(np.float32))
    adata = AnnData(counts.copy())
    adata.var_names = [f"Gene{i}" for i in range(adata.n_vars)]
    adata.obs_names = [f"Cell{i}" for i in range(adata.n_obs)]
    adata.layers["counts"] = counts.copy()
    init_file_registry(adata, "test")

    cfg = {
        "run": {"run_id": "test", "random_seed": 0},
        "keys": {"batch_candidates": []},
        "qc": {"counts_layer": "counts"},
        "core": {
            "engine": "fastcore",
            "fallback_engine": "scanpy_legacy",
            "batch_correction": {"method": "harmony2", "max_iter_harmony": 20},
            "fastcore": {
                "fallback_backend": "scanpy_legacy",
                "allowed_backends": ["fastcore_cpu"],
                "enable_fastcore_cpu_backend": True,
                "omicverse": {
                    "target_sum": 500000,
                    "n_hvgs": 25,
                    "n_pcs": 8,
                    "neighbors": {"n_neighbors": 5},
                    "umap": {"min_dist": 0.3},
                    "leiden": {"resolutions": [0.6], "default_resolution": 0.6},
                },
            },
        },
    }

    result = run_core_pipeline(adata, cfg, tmp_path)

    assert result["backend"] == "fastcore_cpu"
    assert result["fallback_used"] is False
    assert result["batch_correction_method"] == "harmony2"
    assert "pca_covariance_eigh" in result["timings"]
    assert "log1p_norm" in adata.layers
    assert "X_pca_biology" in adata.obsm
    assert "X_pca_identity_prebatch" in adata.obsm
    assert "X_pca_harmony_identity" in adata.obsm
    assert "connectivities_identity" in adata.obsp
    assert "X_umap_identity" in adata.obsm
    assert "cluster_identity" in adata.obs
    assert adata.uns["fastcore_omicverse_gpl"]["license"] == "GPL-3.0-or-later"
    assert (tmp_path / "02_core" / "fastcore" / "fastcore_manifest.json").exists()
