import numpy as np
import pandas as pd
from anndata import AnnData

from fastcore.backends import scanpy_legacy as legacy


def test_legacy_core_reaches_keyword_only_leiden_api(tmp_path, monkeypatch):
    obj = AnnData(np.ones((6, 4)))
    for name in ("normalize_log1p", "score_programs", "robust_z_by_sample", "neighbors_umap"):
        monkeypatch.setattr(legacy, name, lambda *a, **kw: None)
    for name in ("select_hvg", "build_identity_hvg_from_program_decision"):
        monkeypatch.setattr(legacy, name, lambda *a, **kw: pd.DataFrame({"gene": ["g"]}))
    def pca(adata, *, obsm_key, **kwargs):
        adata.obsm[obsm_key] = np.ones((adata.n_obs, 2))
    monkeypatch.setattr(legacy, "run_pca", pca)
    received = []
    def sweep(adata, *, graph_prefix, resolutions, seeds, search_config):
        assert adata is obj
        assert graph_prefix == "identity"
        assert search_config == {"verified_config": True}
        received.append(True)
        return pd.DataFrame({"resolution": resolutions}), pd.DataFrame({"seed": seeds})
    monkeypatch.setattr(legacy, "leiden_sweep", sweep)
    monkeypatch.setattr(legacy, "leiden_search_options", lambda cfg, root: {"verified_config": True})
    result = legacy.run_scanpy_legacy_core(obj, {"core": {"batch_correction": {"method": "none"},
        "leiden_resolutions": [.5], "leiden_seeds": [1, 2]}}, tmp_path)
    assert received == [True] and result["backend"] == "scanpy_legacy"
    assert (tmp_path / "core/tables/leiden_sweep.parquet").is_file()


def test_legacy_core_honours_neighbors_n_pcs(tmp_path, monkeypatch):
    obj = AnnData(np.ones((6, 4)))
    for name in ("normalize_log1p", "score_programs", "robust_z_by_sample"):
        monkeypatch.setattr(legacy, name, lambda *a, **kw: None)
    for name in ("select_hvg", "build_identity_hvg_from_program_decision"):
        monkeypatch.setattr(legacy, name, lambda *a, **kw: pd.DataFrame({"gene": ["g"]}))
    def pca(adata, *, obsm_key, **kwargs):
        adata.obsm[obsm_key] = np.ones((adata.n_obs, 2))
    graphs = {}
    def neighbors_umap(adata, *, prefix, n_pcs, **kwargs):
        graphs[prefix] = n_pcs
    monkeypatch.setattr(legacy, "run_pca", pca)
    monkeypatch.setattr(legacy, "neighbors_umap", neighbors_umap)
    monkeypatch.setattr(legacy, "leiden_sweep", lambda *a, **kw: (pd.DataFrame(), pd.DataFrame()))
    legacy.run_scanpy_legacy_core(obj, {"core": {"batch_correction": {"method": "none"}, "neighbors_n_pcs": 2}}, tmp_path)
    assert graphs == {"identity": 2, "biology": 2}
