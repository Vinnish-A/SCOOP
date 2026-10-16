from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from anndata import AnnData, read_h5ad
from scipy import sparse

from scsp_agent_sop.assays import assay_profile, spatial_qc
from scsp_agent_sop.core import normalize_log1p, select_hvg


def rna():
    counts = np.random.default_rng(21).poisson(np.linspace(0.1, 8, 120), size=(80, 120))
    a = AnnData(sparse.csr_matrix(counts))
    a.var_names = [f"G{i}" for i in range(a.n_vars)]
    a.layers["counts"] = a.X.copy()
    a.obs["sample_id"] = ["A"] * 40 + ["B"] * 40
    normalize_log1p(a)
    return a


def test_rna_default_preserves_seurat_v3_result_and_native_counts(tmp_path):
    a = rna()
    counts = a.layers["counts"].copy()
    expected = sc.pp.highly_variable_genes(
        a, layer="counts", batch_key="sample_id", flavor="seurat_v3", n_top_genes=30, inplace=False
    )
    table = select_hvg(a, n_top_genes=30)
    np.testing.assert_array_equal(a.var.highly_variable_biology, expected.highly_variable)
    assert table.method.unique().tolist() == ["seurat_v3"]
    assert (counts != a.layers["counts"]).nnz == 0
    a.write_h5ad(tmp_path / "rna.h5ad")
    assert not read_h5ad(tmp_path / "rna.h5ad").uns["hvg_selection"]["highly_variable_biology"]["fallback_used"]
    assert assay_profile({}) == "scrna"


def test_loess_fallback_is_clean_log_dispersion_and_records_error(monkeypatch, tmp_path):
    a = rna()
    # Same kind of inherited mixed feature metadata that broke the old batch
    # dispersion fallback. The original metadata should survive unchanged.
    a.var["source_flag"] = pd.Series([True, "native", None] * 40, index=a.var_names, dtype=object)
    original = sc.pp.highly_variable_genes
    calls = []

    def fail_loess(work, **kwargs):
        calls.append((kwargs["flavor"], kwargs.get("batch_key")))
        assert "source_flag" not in work.var
        if kwargs["flavor"] == "seurat_v3":
            raise ValueError(b"There are other near singularities as well.")
        assert kwargs["flavor"] == "seurat"
        assert work.X.max() < 12  # log expression; raw counts are not supplied
        return original(work, **kwargs)

    monkeypatch.setattr(sc.pp, "highly_variable_genes", fail_loess)
    with pytest.warns(RuntimeWarning, match="seurat_dispersion_fallback"):
        table = select_hvg(a, n_top_genes=25)
    assert calls == [("seurat_v3", "sample_id"), ("seurat", "sample_id")]
    assert a.var.source_flag.iloc[1] == "native"
    assert "singularities" in table.fallback_reason.iloc[0]
    info = a.uns["hvg_selection"]["highly_variable_biology"]
    assert info["input_representation"] == "log1p"
    assert json.loads(info["batch_sizes_json"]) == {"A": 40, "B": 40}
    del a.var["source_flag"]  # mixed source metadata itself is not H5AD-serializable
    a.write_h5ad(tmp_path / "fallback.h5ad")


def test_unrelated_hvg_errors_are_not_silently_switched(monkeypatch):
    a = rna()
    monkeypatch.setattr(
        sc.pp, "highly_variable_genes", lambda *a, **k: (_ for _ in ()).throw(ValueError("wrong batch metadata"))
    )
    with pytest.raises(ValueError, match="wrong batch metadata"):
        select_hvg(a)


def test_fractional_v3_counts_are_rejected():
    a = rna()
    a.layers["counts"] = a.layers["counts"].astype(float) / 3
    with pytest.raises(ValueError, match="integer counts"):
        select_hvg(a)


def panel():
    a = AnnData(sparse.csr_matrix([[6, 7, 99], [1, 0, 0], [22, 11, 1], [30, 2, 0]]))
    a.var_names = ["CD1C", "FCER1A", "NegativeControlProbe_1"]
    a.var["feature_types"] = ["Gene Expression", "Gene Expression", "Negative Control Probe"]
    a.layers["counts"] = a.X.copy()
    a.obs["sample_id"] = "A"
    a.obs["cell_area"] = [20.0, 5.0, 0.0, 50.0]
    a.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [np.nan, 3.0]])
    return a


def test_spatial_qc_excludes_controls_and_avoids_rna_rules(tmp_path):
    a = panel()
    cfg = {
        "run": {"assay": "xenium"},
        "qc": {"spatial": {"min_counts": 10, "min_genes": 2, "cell_area_key": "cell_area"}},
    }
    spatial_qc(a, cfg)
    np.testing.assert_array_equal(a.obs.total_counts, [13, 1, 33, 32])
    np.testing.assert_array_equal(a.obs.final_use, [True, False, False, False])
    assert a.obs.pct_counts_mt.isna().all()
    assert a.obs.doublet_call.eq("not_applicable_spatial").all()
    assert a.layers["counts"][0, 2] == 99
    a.write_h5ad(tmp_path / "panel.h5ad")


def test_spatial_gene_only_does_not_claim_zero_background():
    a = panel()[:, :2].copy()
    cfg = {"run": {"assay": "cosmx"}, "qc": {"spatial": {"min_counts": 0, "min_genes": 0}}}
    spatial_qc(a, cfg)
    assert a.obs.non_gene_feature_fraction.isna().all()
    assert not a.uns["qc_assay"]["feature_controls_available"]


def test_panel_requires_feature_semantics():
    a = panel()
    del a.var["feature_types"]
    with pytest.raises(ValueError, match="feature types"):
        spatial_qc(a, {"run": {"assay": "xenium"}})


def test_spatial_core_routes_before_fastcore_and_preserves_full_counts(tmp_path, monkeypatch):
    from scsp_agent_sop.core_runner import run_core_pipeline
    from scsp_agent_sop.storage import init_file_registry

    a = rna()
    a.var["feature_types"] = ["Gene Expression"] * 119 + ["Negative Control Probe"]
    a.obsm["spatial"] = np.random.default_rng(3).normal(size=(80, 2))
    init_file_registry(a, "test")
    before = a.layers["counts"].copy()
    seen = {}

    # Graph construction/PCA run for real; sweep parallelism is independently
    # exercised on actual RA graphs in scripts/benchmark_leiden.py.
    def sweep(adata, **kwargs):
        seen.update(kwargs)
        adata.obs["cluster_identity"] = pd.Categorical(["0"] * adata.n_obs)
        return pd.DataFrame({"resolution": [0.5], "seed": [0]}), pd.DataFrame({"resolution": [0.5], "chosen": [True]})

    monkeypatch.setattr("fastcore.backends.spatial_panel.leiden_sweep", sweep)
    monkeypatch.setattr(
        "scsp_agent_sop.core_runner.plan_fastcore_backend",
        lambda *a, **k: pytest.fail("RNA planner must not run for panel data"),
    )
    cfg = {"run": {"assay": "cosmx"}, "core": {"n_pcs": 6, "spatial": {"min_cells": 1}}}
    result = run_core_pipeline(a, cfg, tmp_path)
    assert result["backend"] == "spatial_scanpy"
    assert not a.var.highly_variable_biology.iloc[-1]
    assert not a.layers["log1p_norm"][:, -1].nnz
    assert (before != a.layers["counts"]).nnz == 0
    assert "stress_score" not in a.obs
    assert seen["search_config"]["backend"] == "networkit_parallel_leiden"
    a.write_h5ad(tmp_path / "panel_core.h5ad")
