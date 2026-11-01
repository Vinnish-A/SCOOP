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


def test_a_configured_symbol_key_missing_from_var_is_an_error():
    import pytest

    obj, _ = _ensembl_indexed()
    with pytest.raises(KeyError, match="qc.symbol_key 'symbol'"):
        annotate_gene_flags(obj, symbol_key="symbol")


def test_qc_keeps_inherited_decisions_as_prior_columns_and_never_readmits(tmp_path, monkeypatch):
    import json
    import sys

    import yaml
    from anndata import read_h5ad
    from scoop.commands import qc

    rng = np.random.default_rng(32)
    counts = sparse.csr_matrix(rng.poisson(rng.gamma(1.0, 2.0, size=(200, 120))).astype(np.float32))
    obj = ad.AnnData(counts.copy(), obs=pd.DataFrame({"lib": ["a"] * 100 + ["b"] * 100},
                                                     index=[f"c{i}" for i in range(200)]))
    obj.var_names = [f"G{i}" for i in range(118)] + ["MT-ND1", "MT-CO1"]
    obj.layers["counts"] = counts
    obj.obs["final_use"] = np.arange(200) >= 30  # an earlier ambient/sample mask excluded 30 cells
    obj.obs["qc_class"] = "fail"
    obj.obs["total_counts"] = 99999.0
    run = tmp_path / "runs" / "r1"
    (run / "config").mkdir(parents=True)
    (run / "artifacts").mkdir()
    obj.write_h5ad(run / "artifacts" / "validated_adata.h5ad")
    config = run / "config" / "run.yaml"
    config.write_text(yaml.safe_dump({"run": {"obs_keys": {"sample": "lib"}}, "qc": {"scrublet": {"workers": 1}}}))
    monkeypatch.setattr(sys, "argv", ["scoop qc", "--config", str(config)])
    qc.main()

    out = read_h5ad(run / "artifacts" / "adata_qc.h5ad")
    assert not out.obs.final_use.iloc[:30].any()
    assert (out.obs.final_use <= out.obs.prior_final_use).all()
    assert out.obs.prior_qc_class.eq("fail").all() and out.obs.prior_total_counts.eq(99999.0).all()
    np.testing.assert_allclose(out.obs.total_counts, np.asarray(counts.sum(axis=1)).ravel())
    decision = json.loads((run / "logs" / "decision_log.jsonl").read_text().splitlines()[-1])
    assert decision["parameters"]["sample_key"] == "lib"
    assert decision["evidence"]["prior_qc_columns"] == {"final_use": "prior_final_use", "qc_class": "prior_qc_class",
                                                        "total_counts": "prior_total_counts"}
    assert decision["evidence"]["prior_excluded_cells"] == 30
