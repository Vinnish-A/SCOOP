import json

import numpy as np
import pandas as pd
import pytest
import yaml
from anndata import AnnData, read_h5ad
from scipy import sparse

from scoop.atlas import _match_index_dtypes, annotate, integrate


def _run(tmp_path, scvi=None, split=None, extra=None, graph=None):
    rng = np.random.default_rng(0)
    counts = sparse.csr_matrix(rng.poisson(2.0, (60, 20)).astype(np.float32))
    obs = pd.DataFrame({"donor_id": np.repeat(["d1", "d2", "d3"], 20), "batch": np.tile(["b1", "b2"], 30)},
                       index=[f"c{i}" for i in range(60)])
    var = pd.DataFrame({"hv": np.arange(20) < 12}, index=[f"g{i}" for i in range(20)])
    source = tmp_path / "in.h5ad"
    AnnData(counts, obs=obs, var=var, layers={"counts": counts}).write_h5ad(source)
    root = tmp_path / "runs/r1"
    (root / "config").mkdir(parents=True, exist_ok=True)
    config = root / "config/run.yaml"
    config.write_text(yaml.safe_dump({"atlas": {
        "split": {"unit": "donor_id", **(split or {})},
        "scvi": {"n_latent": 2, "max_epochs": 1, "batch_size": 32, "accelerator": "cpu", "devices": 1, **(scvi or {})},
        "graph": {"n_neighbors": 5, "leiden_resolutions": [1.0], "leiden_seeds": [0],
                  "leiden_search": {"backend": "scanpy"}, **(graph or {})}}, **(extra or {})}))
    return integrate(config_path=config, input_path=source), root


def test_scvi_batch_key_must_be_declared(tmp_path):
    with pytest.raises(ValueError, match="atlas.scvi.batch_key"):
        _run(tmp_path, scvi={"batch_key": None})


def test_held_out_units_are_projected_not_trained(tmp_path):
    out, root = _run(tmp_path, scvi={"batch_key": "batch"}, split={"held_out_query_units": ["d3"]})
    adata = read_h5ad(out)
    assert adata.obs.is_projected.tolist() == [False] * 40 + [True] * 20
    assert adata.obsm["X_scVI"].shape == (60, 2) and np.isfinite(adata.obsm["X_scVI"]).all()
    params = json.loads((root / "atlas/scvi-run.json").read_text())["parameters"]
    assert params["batch_key"] == "batch" and params["held_out_query_units"] == ["d3"]
    assert (params["n_training_cells"], params["n_projected_cells"]) == (40, 20)


def test_held_out_units_must_exist_and_share_training_batches(tmp_path):
    with pytest.raises(ValueError, match="not found"):
        _run(tmp_path, scvi={"batch_key": "batch"}, split={"held_out_query_units": ["d9"]})
    with pytest.raises(ValueError, match="absent from the training cells"):
        _run(tmp_path, scvi={"batch_key": "donor_id"}, split={"held_out_query_units": ["d3"]})


def test_typed_batch_key_is_accepted_but_a_guessed_candidate_is_not(tmp_path):
    (tmp_path / "typed").mkdir()
    (tmp_path / "guessed").mkdir()
    result, root = _run(tmp_path / "typed", scvi={"batch_key": None}, extra={"run": {"obs_keys": {"batch": "batch"}}})
    assert '"batch_key_source": "run.obs_keys.batch"' in (root / "atlas/scvi-run.json").read_text()
    with pytest.raises(ValueError, match="run.obs_keys.batch"):
        _run(tmp_path / "guessed", scvi={"batch_key": None}, extra={"keys": {"batch_candidates": ["batch"]}})


def test_declared_covariate_must_exist(tmp_path):
    with pytest.raises(KeyError, match="chemistry"):
        _run(tmp_path, scvi={"batch_key": "batch", "categorical_covariate_keys": ["chemistry"]})


def test_model_and_history_survive_a_late_graph_failure(tmp_path, monkeypatch):
    import scanpy as sc

    def fail(*args, **kwargs):
        raise RuntimeError("simulated UMAP failure")

    monkeypatch.setattr(sc.tl, "umap", fail)
    with pytest.raises(RuntimeError, match="simulated UMAP"):
        _run(tmp_path, scvi={"batch_key": "batch"})
    root = tmp_path / "runs/r1"
    assert (root / "atlas/models/scvi/model.pt").is_file()
    assert list((root / "atlas/tables").glob("training_*.tsv"))


def test_model_kwargs_hvg_key_and_chosen_resolution_are_recorded(tmp_path):
    out, root = _run(tmp_path, scvi={"batch_key": "batch", "hvg_key": "hv", "n_latent": 3,
                                     "model_kwargs": {"n_layers": 2, "n_hidden": 16, "dispersion": "gene-batch"}},
                     split={"held_out_query_units": ["d3"]}, graph={"leiden_resolutions": [0.2, 2.0]})
    adata = read_h5ad(out)
    params = json.loads((root / "atlas/scvi-run.json").read_text())["parameters"]

    assert params["model_kwargs"] == {"n_layers": 2, "n_hidden": 16, "dispersion": "gene-batch", "n_latent": 3,
                                      "gene_likelihood": "nb"}
    assert (params["hvg_key"], params["n_training_genes"], adata.n_vars) == ("hv", 12, 20)
    assert adata.obsm["X_scVI"].shape == (60, 3)
    stability = pd.read_csv(root / "atlas/tables/leiden_stability.tsv", sep="\t")
    assert params["leiden_resolution"] == float(stability.loc[stability.chosen, "resolution"].iloc[0])
    assert params["leiden_resolution"] == adata.uns["leiden_selection"]["leiden_scvi"]["chosen_resolution"]


def test_mixed_int64_indptr_is_cast_to_int32():
    matrix = sparse.csr_matrix(np.eye(4, dtype=np.float32))
    matrix.indptr = matrix.indptr.astype(np.int64)
    adata = AnnData(matrix.copy(), layers={"counts": matrix})
    adata.X.indptr = adata.X.indptr.astype(np.int64)

    _match_index_dtypes(adata)

    for value in (adata.X, adata.layers["counts"]):
        assert (value.indptr.dtype, value.indices.dtype) == (np.int32, np.int32)
        np.testing.assert_array_equal(value.toarray(), np.eye(4))


def test_annotate_tests_log1p_cp10k_and_registers_a_label_keyed_deg_table(tmp_path):
    out, root = _run(tmp_path, scvi={"batch_key": "batch"})

    target, _, evidence_path = annotate(config_path=root / "config/run.yaml", input_path=out)

    adata = read_h5ad(target)
    counts = adata.layers["counts"].toarray()
    expected = np.log1p(counts / counts.sum(axis=1, keepdims=True) * 1e4)
    np.testing.assert_allclose(adata.layers["log1p_norm"].toarray(), expected, rtol=1e-5)
    np.testing.assert_array_equal(adata.X.toarray(), counts)
    record = adata.uns["file_registry"]["tables"]["scvi_cluster_markers"]
    assert record["label_key"] == "leiden_scvi"
    assert not [key for key in adata.uns if key.startswith("top_markers_")]
    clusters = json.loads(evidence_path.read_text())["clusters"]
    assert clusters and all(c["deg_refs"] == ["scvi_cluster_markers"] and c["top_markers"] for c in clusters)
