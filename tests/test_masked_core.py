import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import sparse

from scsp_agent_sop.core import _masked_hvg, normalize_log1p, score_programs, select_hvg
from scsp_agent_sop.measurement import FeatureMeasurement


def _rna(n_genes=120, seed=21):
    rng = np.random.default_rng(seed)
    counts = rng.poisson(np.linspace(0.3, 8, n_genes), size=(120, n_genes)).astype(float)
    a = AnnData(sparse.csr_matrix(counts))
    a.var_names = [f"G{i}" for i in range(n_genes)]
    a.obs["sample_id"] = ["A"] * 40 + ["B"] * 40 + ["C"] * 40
    a.layers["counts"] = a.X.copy()
    return a


def _unmeasured_in_b(a, gene):
    """``gene`` was not quantified in sample B: its counts there are structural zeros."""
    j = a.var_names.get_loc(gene)
    counts = a.layers["counts"].toarray()
    counts[40:80, j] = 0
    a.layers["counts"] = sparse.csr_matrix(counts)
    measured = np.ones((3, a.n_vars), bool)
    measured[1, j] = False
    return FeatureMeasurement.from_var("sample_id", ["A", "B", "C"], a.var_names, measured)


def test_masked_merge_reproduces_scanpy_batches_when_everything_is_measured():
    a = _rna()
    batches = a.obs.sample_id.to_numpy()
    for flavor in ("seurat_v3", "seurat_v3_paper"):
        expected = sc.pp.highly_variable_genes(a, flavor=flavor, batch_key="sample_id", n_top_genes=30, inplace=False)
        got = _masked_hvg(a, batches, pd.unique(batches), np.ones((3, a.n_vars), bool), flavor=flavor, n_top_genes=30)
        np.testing.assert_array_equal(got.highly_variable, expected.highly_variable)
        for col in ("means", "variances", "variances_norm", "highly_variable_rank"):
            np.testing.assert_allclose(got[col], expected[col], rtol=1e-6, equal_nan=True)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    expected = sc.pp.highly_variable_genes(a, flavor="seurat", batch_key="sample_id", n_top_genes=30, inplace=False)
    got = _masked_hvg(a, batches, pd.unique(batches), np.ones((3, a.n_vars), bool), flavor="seurat", n_top_genes=30)
    np.testing.assert_array_equal(got.highly_variable, expected.highly_variable)
    np.testing.assert_array_equal(got.highly_variable_intersection, expected.highly_variable_intersection)


def test_fully_measured_measurement_keeps_the_plain_result():
    plain, masked = _rna(), _rna()
    fm = FeatureMeasurement.from_var("sample_id", ["A", "B", "C"], masked.var_names, np.ones((3, masked.n_vars), bool))
    expected = select_hvg(plain, n_top_genes=30)
    got = select_hvg(masked, n_top_genes=30, measurement=fm)
    pd.testing.assert_frame_equal(got, expected)
    assert "measurement_sample_key" not in plain.uns["hvg_selection"]["highly_variable_biology"]
    assert masked.uns["hvg_selection"]["highly_variable_biology"]["n_genes_partially_measured"] == 0


def test_unmeasured_gene_has_batch_statistics_from_measuring_cells_only():
    a = _rna()
    fm = _unmeasured_in_b(a, "G100")
    table = select_hvg(a, n_top_genes=30, measurement=fm).set_index("gene")
    measuring = a.layers["counts"][np.r_[0:40, 80:120], a.var_names.get_loc("G100")].toarray().ravel()
    assert table.loc["G100", "n_batches_measured"] == 2
    np.testing.assert_allclose(table.loc["G100", "means"], measuring.mean())
    np.testing.assert_allclose(table.loc["G100", "variances"], measuring.var(ddof=1))
    assert table.loc["G5", "n_batches_measured"] == 3
    assert a.uns["hvg_selection"]["highly_variable_biology"]["n_genes_partially_measured"] == 1


def test_single_batch_never_ranks_a_gene_on_structural_zeros():
    a = _rna()
    j = a.var_names.get_loc("G3")
    counts = a.layers["counts"].toarray()
    counts[:, j] = np.r_[np.full(40, 10.0), np.zeros(40), np.full(40, 10.0)]  # bimodal only through B's zeros
    a.layers["counts"] = sparse.csr_matrix(counts)
    fm = _unmeasured_in_b(a, "G3")
    plain = select_hvg(a.copy(), batch_key=None, n_top_genes=10)
    masked = select_hvg(a, batch_key=None, n_top_genes=10, measurement=fm).set_index("gene")
    assert plain.set_index("gene").loc["G3", "highly_variable"]
    assert not masked.loc["G3", "highly_variable"] and masked.loc["G3", "n_batches_measured"] == 0
    assert masked.highly_variable.sum() == 10


def _scoring_object():
    stress = ["FOS", "JUN", "JUNB", "JUND", "ATF3", "EGR1", "DUSP1", "HSPA1A", "HSPA1B", "DNAJB1"]
    cycle = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TOP2A", "MKI67"]
    names = stress + cycle + ["RPL3", "RPS6"] + [f"F{i}" for i in range(276)]
    a = _rna(n_genes=len(names), seed=3)
    a.var_names = names
    a.var["ribo_gene"] = a.var_names.str.startswith(("RPL", "RPS"))
    return a


def test_scores_ignore_genes_a_sample_did_not_measure():
    a, b = _scoring_object(), _scoring_object()
    fm = _unmeasured_in_b(a, "FOS")
    counts = a.layers["counts"].toarray()
    counts[:, a.var_names.get_loc("FOS")] = 0.0  # whatever FOS holds never enters any score
    b.layers["counts"] = sparse.csr_matrix(counts)
    # Plain log1p: library-size normalisation would let FOS shift every other gene.
    for obj in (a, b):
        obj.layers["log1p_norm"] = sparse.csr_matrix(np.log1p(obj.layers["counts"].toarray()))
        score_programs(obj, measurement=fm)
    columns = ["stress_score", "cell_cycle_s_score", "cell_cycle_g2m_score", "ribo_score"]
    pd.testing.assert_frame_equal(a.obs[columns], b.obs[columns])
    assert a.uns["program_score_genes"]["stress_score"] == {"n_on_axis": 10, "n_used": 9}
    assert a.uns["program_score_genes"]["ribo_score"] == {"n_on_axis": 2, "n_used": 2}


def test_scores_without_measurement_are_scanpys():
    a = _scoring_object()
    normalize_log1p(a)
    expected = a.copy()
    score_programs(a)
    expected.X = expected.layers["log1p_norm"].copy()
    sc.tl.score_genes(expected, list(a.var_names[:10]), score_name="stress_score")
    np.testing.assert_array_equal(a.obs.stress_score, expected.obs.stress_score)
    assert "program_score_genes" not in a.uns


def test_rescoring_replaces_cell_cycle_columns():
    a = _rna()
    a.var_names = ["FOS", "JUN", "ATF3", "MCM5", "PCNA", "TYMS", "TOP2A", "MKI67", "CDK1"] + list(a.var_names[9:])
    normalize_log1p(a)
    score_programs(a)
    score_programs(a)
    assert list(a.obs.columns).count("cell_cycle_s_score") == 1
    assert list(a.obs.columns).count("cell_cycle_g2m_score") == 1


def test_size_factor_key_normalises_a_gene_subset_by_full_library_size():
    import pytest

    full = _rna()
    full.obs["total_counts"] = np.asarray(full.layers["counts"].sum(axis=1)).ravel()
    subset = full[:, :30].copy()
    normalize_log1p(subset, size_factor_key="total_counts")
    counts = subset.layers["counts"].toarray()
    expected = np.log1p(counts / full.obs.total_counts.to_numpy()[:, None] * 1e4)
    np.testing.assert_allclose(subset.layers["log1p_norm"].toarray(), expected, rtol=1e-6)
    normalize_log1p(full, size_factor_key="total_counts")
    plain = _rna()
    normalize_log1p(plain)
    np.testing.assert_allclose(full.layers["log1p_norm"].toarray(), plain.layers["log1p_norm"].toarray(), rtol=1e-6)
    subset.obs["log10_total"] = np.log10(subset.obs.total_counts)
    with pytest.raises(ValueError, match="no smaller than"):
        normalize_log1p(subset, size_factor_key="log10_total")


def test_a_gene_set_absent_from_the_axis_is_unassessed_not_a_crash_or_zero():
    from scsp_agent_sop.core import robust_z_by_sample

    a = _rna()  # a panel without stress, cell-cycle or ribosomal genes
    normalize_log1p(a)
    score_programs(a)
    for column in ("stress_score", "ribo_score", "cell_cycle_s_score", "cell_cycle_g2m_score", "proliferation_score"):
        assert a.obs[column].isna().all()
    assert a.uns["program_score_genes"]["stress_score"] == {"n_on_axis": 0, "n_used": 0}
    robust_z_by_sample(a, ["stress_score"], sample_key="sample_id")
    assert a.obs.stress_score_z_sample.isna().all()
