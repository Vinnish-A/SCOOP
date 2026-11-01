from types import SimpleNamespace

import numpy as np
import pytest
import json
import pandas as pd
from anndata import AnnData, read_h5ad
from scipy import sparse

from scsp_agent_sop.rctd_adapter import selected_weights
from scsp_agent_sop import rctd_adapter as adapter
from scoop.io import sha256_file


def test_multi_uses_sparse_solution_not_initial_full_fit():
    result = SimpleNamespace(
        weights=np.full((2, 3), 0.3),
        sub_weights=np.array([[0.8, 0.2], [1.0, 0.0]]),
        cell_type_indices=np.array([[2, 0], [1, -1]]),
    )
    np.testing.assert_allclose(selected_weights(result, "multi"), [[0.2, 0.0, 0.8], [0.0, 1.0, 0.0]])
    np.testing.assert_array_equal(selected_weights(result, "full"), result.weights)


def test_doublet_preserves_upstream_pair_weights_separately_from_diagnostic_class():
    result = SimpleNamespace(
        weights=np.full((2, 3), 0.3),
        first_type=np.array([0, 1]),
        second_type=np.array([2, 1]),
        weights_doublet=np.array([[0.7, 0.3], [0.6, 0.4]]),
    )
    np.testing.assert_allclose(selected_weights(result, "doublet"), [[0.7, 0.0, 0.3], [0.0, 1.0, 0.0]])
    with pytest.raises(ValueError):
        selected_weights(result, "unknown")


def test_source_exclusions_and_native_exposure_precede_gene_projection():
    obj = AnnData(sparse.csr_matrix([[30, 70, 900], [40, 60, 0], [25, 75, 0]]),
        obs=pd.DataFrame({"source_author_retained": pd.array([True, False, pd.NA], dtype="boolean")}, index=["a", "b", "c"]),
        var=pd.DataFrame({"is_gene_expression": [True, True, False]}, index=["g1", "g2", "control"]))
    full, status = adapter._count_input(obj, {}, {"spatial_counts_layer": "X"})
    assert full.shape == (3, 2)
    assert status.source_selected.tolist() == [True, False, True]
    assert status.source_selection_partial.tolist() == [False, False, True]
    assert status.native_full_exposure.tolist() == [100, 100, 100]


def test_count_input_rejects_normalized_expression_and_registry_axis_mismatch(tmp_path):
    obj = AnnData(np.array([[1.2, 0.8]]), obs=pd.DataFrame(index=["a"]), var=pd.DataFrame(index=["g1", "g2"]))
    with pytest.raises(ValueError, match="integer counts"):
        adapter._count_input(obj, {}, {"spatial_counts_layer": "X"})
    registry = tmp_path / "axis.parquet"
    pd.DataFrame({"source_axis_id": ["different"]}).to_parquet(registry)
    with pytest.raises(ValueError, match="IDs/order"):
        adapter._count_input(obj, {}, {"spatial_counts_layer": "X", "spatial_registry": str(registry)})


def test_reference_retains_rare_type_and_full_depth_when_excluding_genes(tmp_path):
    obj = AnnData(sparse.csr_matrix([[70, 30, 100], [20, 80, 100], [10, 90, 100]]),
        obs=pd.DataFrame({"cell_type_lvl3": ["rare", "common", "common"], "annotation_status": ["accepted"] * 3}, index=["a", "b", "c"]),
        var=pd.DataFrame(index=["heldout", "g2", "g3"]))
    path = tmp_path / "reference.h5ad"; obj.write_h5ad(path)
    ref, info = adapter._reference(path, {}, {"reference_counts_layer": "X", "exclude_genes": ["heldout"]})
    assert set(ref.cell_type_names) == {"rare", "common"}
    assert info["rare_components"] == ["common", "rare"]
    np.testing.assert_allclose(ref.profiles[:, ref.cell_type_names.index("rare")], [.15, .5])
    with pytest.raises(ValueError, match="fewer than|dropped declared"):
        adapter._reference(path, {}, {"reference_counts_layer": "X", "reference_cell_min": 2})


def test_projected_reference_cannot_silently_change_native_denominator(tmp_path):
    obj = AnnData(np.array([[60, 40], [40, 60]]),
        obs=pd.DataFrame({"cell_type_lvl3": ["a", "b"], "annotation_status": ["accepted"] * 2,
                          "full_depth": [1000, 1000]}, index=["a", "b"]), var=pd.DataFrame(index=["g1", "g2"]))
    path = tmp_path / "reference.h5ad"; obj.write_h5ad(path)
    with pytest.raises(ValueError, match="Projected reference"):
        adapter._reference(path, {}, {"reference_counts_layer": "X", "reference_exposure_key": "full_depth"})


def test_valid_fit_survives_export_failure_without_refitting(tmp_path, monkeypatch):
    obj = AnnData(sparse.csr_matrix([[70, 30], [20, 80], [30, 70]]),
        obs=pd.DataFrame({"source_author_retained": [True, False, True]}, index=["a", "b", "c"]),
        var=pd.DataFrame(index=["g1", "g2"]))
    path = tmp_path / "spatial.h5ad"; obj.write_h5ad(path)
    profiles = tmp_path / "profiles.npz"
    np.savez(profiles, profiles=np.array([[.7, .2], [.3, .8]]), gene_names=np.array(["g1", "g2"]), cell_type_names=np.array(["A", "B"]))
    manifest = tmp_path / "reference.json"
    manifest.write_text(json.dumps({"schema_version": "scoop.reference-profiles.v1", "unit": "mean_counts_per_full_native_exposure",
        "labels_reviewed": True, "native_full_exposure_preserved": True, "profiles_sha256": sha256_file(profiles),
        "training_cells": {"A": 20, "B": 20}}))
    lookup = tmp_path / "q.npz"; lookup.write_bytes(b"test fixture lookup")
    cfg = {"spatial": {"rctd": {"spatial_counts_layer": "X", "reference_profiles": str(profiles),
        "reference_manifest": str(manifest), "q_matrices": str(lookup)}}}
    calls = []
    def fake_fit(spatial, status, reference, options, mode):
        calls.append(True)
        assert status.source_selected.tolist() == [True, False, True]
        return {"raw": np.array([[.8, .2], [.3, .7]]), "fitted_rows": np.array([0, 2]),
                "profiles": reference.profiles, "fit_genes": np.array(reference.gene_names),
                "types": np.array(reference.cell_type_names), "converged": np.array([True, True])}, {
                "method": "test fixture", "version": "0.4.0", "mode": mode, "sigma": 100}
    monkeypatch.setattr(adapter, "_fit_payload", fake_fit)
    writer = adapter.write_h5ad_atomic
    def fail_export(*args, **kwargs):
        raise OSError("serialization failed after fit")
    monkeypatch.setattr(adapter, "write_h5ad_atomic", fail_export)
    output = tmp_path / "out/result.h5ad"
    with pytest.raises(OSError, match="serialization"):
        adapter.run_rctd_api(path, None, output, cfg, "full")
    assert len(list((tmp_path / "out/cache").rglob("fit.json"))) == 1
    monkeypatch.setattr(adapter, "write_h5ad_atomic", writer)
    info = adapter.run_rctd_api(path, None, output, cfg, "full")
    assert len(calls) == 1 and info["computation_reused"]
    saved = read_h5ad(output)
    assert saved.obs_names.tolist() == ["a", "b", "c"]
    assert np.isnan(saved.obsm["rctd_weights_raw"][1]).all()
    np.testing.assert_allclose(saved.obsm["rctd_weights_raw"][[0, 2]], [[.8, .2], [.3, .7]])
