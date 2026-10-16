from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd

from scsp_agent_sop.storage import prune_h5ad


def test_prune_h5ad_minimal_schema():
    a = ad.AnnData(np.ones((2, 2)))
    a.obs["sample_id"] = ["s1", "s1"]
    a.obs["tmp"] = [1, 2]
    a.var["gene_symbol"] = ["A", "B"]
    a.var["tmp"] = [1, 2]
    a.uns["file_registry"] = {}
    a.uns["big_tmp"] = {"x": list(range(10))}
    schema = {"obs_keep": {"identity": ["sample_id"]}, "var_keep": ["gene_symbol"], "layers_keep": [], "obsm_keep": [], "obsp_keep": [], "uns_keep": ["file_registry"]}
    prune_h5ad(a, schema)
    assert "sample_id" in a.obs
    assert "tmp" not in a.obs
    assert "gene_symbol" in a.var
    assert "big_tmp" not in a.uns


def test_delivery_retains_executable_graph_descriptors_and_umap(tmp_path):
    from scipy import sparse
    obj = ad.AnnData(np.ones((4, 3)))
    obj.obsp['connectivities_identity'] = sparse.eye(4, format='csr')
    obj.obsp['distances_identity'] = sparse.eye(4, format='csr')
    obj.uns['neighbors_identity'] = {'connectivities_key': 'neighbors_identity_connectivities',
                                     'distances_key': 'neighbors_identity_distances', 'params': {'n_neighbors': 3}}
    obj.uns['umap'] = {'params': {'a': 1., 'b': 1.}}
    obj.obsm['X_umap_identity'] = np.ones((4, 2))
    obj.uns['irrelevant'] = {'large': np.ones(50)}
    prune_h5ad(obj, {'obsp_keep': ['connectivities_identity', 'distances_identity'],
                     'obsm_keep': ['X_umap_identity'], 'uns_keep': []})
    path = tmp_path / 'delivery.h5ad'
    obj.write_h5ad(path)
    reloaded = ad.read_h5ad(path)
    assert reloaded.uns['neighbors_identity']['connectivities_key'] == 'connectivities_identity'
    assert reloaded.uns['neighbors_identity']['distances_key'] == 'distances_identity'
    assert 'umap' in reloaded.uns
    assert 'irrelevant' not in reloaded.uns


def test_atomic_export_preserves_previous_object_on_serialization_failure(tmp_path):
    import pytest
    from scsp_agent_sop.storage import write_h5ad_atomic, register_file
    obj = ad.AnnData(np.ones((2, 2)))
    path = tmp_path / 'working.h5ad'
    write_h5ad_atomic(obj, path)
    before = path.read_bytes()
    obj.uns['invalid_records'] = [{'nested': 1}, {'nested': 2}]
    with pytest.raises((TypeError, ValueError)):
        write_h5ad_atomic(obj, path)
    assert path.read_bytes() == before
    del obj.uns['invalid_records']
    register_file(obj, key='self', path=path, category='artifacts', fingerprint=False)
    write_h5ad_atomic(obj, path)
    assert 'sha256' not in ad.read_h5ad(path).uns['file_registry']['artifacts']['self']
