from pathlib import Path
import anndata as ad
import numpy as np
import pytest

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


def test_prune_h5ad_keeps_only_what_the_schema_names():
    a = ad.AnnData(np.ones((2, 2)))
    a.obs["annotation_reviewed"] = [True, True]
    a.obs["cell_type_lvl3"] = ["T", "B"]
    a.uns["annotation_review"] = {"reviewed": 2}
    a.uns["delivery_scope"] = {"scope": "stale"}
    prune_h5ad(a, {"obs_keep": {"clustering_annotation": ["cell_type_lvl3"]}, "uns_keep": ["annotation_review"]})
    assert list(a.obs.columns) == ["cell_type_lvl3"] and list(a.uns) == ["annotation_review"]


@pytest.mark.parametrize("schema", ["configs/h5ad_schema.yaml", "src/scoop/data/h5ad-schema.yaml"])
def test_delivery_schema_retains_review_cnv_projection_and_programme_fields(schema):
    import yaml
    a = ad.AnnData(np.ones((2, 2)))
    obs = ["annotation_reviewed", "is_projected", "fastcnv_reference_pool", "fastcnv_cnv_fraction",
           "fastcnv_normal_threshold", "fastcnv_threshold_source", "fastcnv_tumor_evidence"]
    for column in obs:
        a.obs[column] = ["x", "y"]
    a.var["nmf_discovery_gene"] = [True, False]
    a.uns["annotation_review"] = {"reviewed": 2}
    prune_h5ad(a, yaml.safe_load((Path(__file__).resolve().parents[1] / schema).read_text()))
    assert list(a.obs.columns) == obs and list(a.var.columns) == ["nmf_discovery_gene"] and "annotation_review" in a.uns


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


def test_atomic_writers_apply_the_process_umask(tmp_path):
    import os
    from scoop.io import atomic_json, atomic_text
    from scsp_agent_sop.config import resolve_run_root
    from scsp_agent_sop.storage import write_h5ad_atomic
    previous = os.umask(0o022)
    try:
        h5ad = write_h5ad_atomic(ad.AnnData(np.ones((2, 2))), tmp_path / 'out.h5ad')
        receipt = atomic_json(tmp_path / 'receipt.json', {'ok': True})
        gate = atomic_text(tmp_path / 'gate.yaml', 'x', mode=0o600)
    finally:
        os.umask(previous)
    assert h5ad.stat().st_mode & 0o777 == 0o644
    assert receipt.stat().st_mode & 0o777 == 0o644
    assert gate.stat().st_mode & 0o777 == 0o600
    # A config outside runs/<id>/config/ roots the run beside it, not under the cwd.
    assert resolve_run_root(tmp_path / 'elsewhere' / 'run.yaml', {'run': {'run_id': 'x'}}) == (tmp_path / 'elsewhere').resolve()
    from scsp_agent_sop.config import project_root_from_run_root
    assert project_root_from_run_root(tmp_path / 'elsewhere') == tmp_path / 'elsewhere'
    assert project_root_from_run_root(tmp_path / 'runs' / 'r1') == tmp_path


def test_default_file_mode_reads_the_umask_without_setting_it(monkeypatch):
    import os
    import sys
    from scoop.io import default_file_mode
    if not sys.platform.startswith('linux'):
        pytest.skip('procfs umask is Linux-only')
    previous = os.umask(0o027)
    try:
        def fail(_):
            raise AssertionError('os.umask changes the mask of every thread')
        monkeypatch.setattr(os, 'umask', fail)
        assert default_file_mode() == 0o640
    finally:
        monkeypatch.undo()
        os.umask(previous)
