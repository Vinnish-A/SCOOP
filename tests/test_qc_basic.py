import warnings

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from scsp_agent_sop.qc import annotate_gene_flags, assign_qc_flags_per_sample, compute_basic_qc


def _ensembl_indexed():
    counts = np.random.default_rng(3).poisson(3, size=(40, 4))
    var = pd.DataFrame({"gene_symbol": ["MT-CO1", "RPL3", "ACTB", "GAPDH"]}, index=[f"ENSG{i:011d}" for i in range(4)])
    obj = ad.AnnData(sparse.csr_matrix(counts), var=var)
    obj.layers["counts"] = obj.X.copy()
    obj.obs["sample_id"] = ["a"] * 20 + ["b"] * 20
    return obj, counts


def test_unflagged_mitochondria_are_unassessed_not_zero():
    obj, _ = _ensembl_indexed()
    annotate_gene_flags(obj)
    compute_basic_qc(obj)
    assert obj.obs.pct_counts_mt.isna().all() and obj.obs.mt_unassessed.all()
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        table = assign_qc_flags_per_sample(obj)
    assert not obj.obs.qc_high_mt_flag.any()
    assert table.pct_counts_mt_median.isna().all()


def test_symbol_key_flags_genes_on_an_ensembl_axis():
    obj, counts = _ensembl_indexed()
    annotate_gene_flags(obj, symbol_key="gene_symbol")
    compute_basic_qc(obj)
    assert not obj.obs.mt_unassessed.any()
    np.testing.assert_allclose(obj.obs.pct_counts_mt, 100 * counts[:, 0] / np.maximum(counts.sum(1), 1))
    assert obj.obs.hb_unassessed.all() and obj.obs.pct_counts_hb.isna().all()


def test_samples_that_measured_no_mt_gene_are_unassessed():
    from scsp_agent_sop.measurement import FeatureMeasurement

    obj, counts = _ensembl_indexed()
    annotate_gene_flags(obj, symbol_key="gene_symbol")
    measured = np.array([[True, True, True, True], [False, True, True, True]])
    measurement = FeatureMeasurement.from_var("sample_id", ["a", "b"], obj.var_names, measured)
    measurement.check(obj)
    compute_basic_qc(obj, measurement=measurement)
    in_b = obj.obs.sample_id.eq("b").to_numpy()
    np.testing.assert_array_equal(obj.obs.mt_unassessed, in_b)
    assert obj.obs.pct_counts_mt[in_b].isna().all()
    np.testing.assert_allclose(obj.obs.pct_counts_mt[~in_b], 100 * counts[~in_b, 0] / np.maximum(counts[~in_b].sum(1), 1))
    assert not obj.obs.ribo_unassessed.any()
