import sys

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import yaml
from anndata import AnnData, read_h5ad
from scipy.stats import mannwhitneyu

from scsp_agent_sop.annotation import _mwu_rows, run_markers_scanpy
from scsp_agent_sop.measurement import FORMAT, FeatureMeasurement
from scsp_agent_sop.storage import register_file


def _adata(n=150, n_genes=6, seed=0):
    rng = np.random.default_rng(seed)
    counts = rng.poisson(rng.uniform(0.3, 3.0, n_genes), (n, n_genes)).astype(float)
    obs = pd.DataFrame({"cluster": pd.Categorical(rng.choice(["a", "b", "c"], n)),
                        "sample": rng.choice(["s1", "s2", "s3"], n)}, index=[f"c{i}" for i in range(n)])
    adata = AnnData(sp.csr_matrix(np.log1p(counts)), obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]))
    measured = np.ones((3, n_genes), bool)
    measured[0, 4] = False  # g4 unmeasured in s1
    measured[[0, 1], 5] = False  # g5 measured in s3 only
    unmeasured = np.isin(obs["sample"], ["s1"])[:, None] & (np.arange(n_genes) == 4) | \
        np.isin(obs["sample"], ["s1", "s2"])[:, None] & (np.arange(n_genes) == 5)
    adata.X = sp.csr_matrix(np.where(unmeasured, 0.0, adata.X.toarray()))  # structural zeros
    return adata, FeatureMeasurement.from_var("sample", ["s1", "s2", "s3"], adata.var_names, measured)


def test_fully_measured_genes_keep_the_scanpy_result():
    adata, fm = _adata()
    masked = run_markers_scanpy(adata, groupby="cluster", layer="X", measurement=fm)
    plain = run_markers_scanpy(adata[:, :4].copy(), groupby="cluster", layer="X")
    merged = plain.merge(masked, on=["group", "names"], suffixes=("_plain", "_masked"))
    assert len(merged) == len(plain) == 12
    for column in ("scores", "pvals", "logfoldchanges", "pct_nz_group", "U", "log10_pval"):
        np.testing.assert_allclose(merged[f"{column}_masked"], merged[f"{column}_plain"], rtol=1e-12)
    # the masked engine run on every cell reproduces scanpy's tie-corrected z
    engine = _mwu_rows(adata.X, None, np.arange(4), adata.obs.cluster.cat.codes.to_numpy(),
                       np.asarray(adata.obs.cluster.cat.categories), adata.var_names, 0)
    check = engine.merge(plain, on=["group", "names"], suffixes=("_engine", ""))
    np.testing.assert_allclose(check.scores_engine, check.scores, atol=1e-5)
    np.testing.assert_allclose(check.logfoldchanges_engine, check.logfoldchanges, atol=1e-5)


def test_partially_measured_genes_use_measuring_cells_only():
    adata, fm = _adata()
    result = run_markers_scanpy(adata, groupby="cluster", layer="X", measurement=fm)
    x = adata.X.toarray()
    for gene, samples in (("g0", ["s1", "s2", "s3"]), ("g4", ["s2", "s3"]), ("g5", ["s3"])):
        cells = adata.obs["sample"].isin(samples).to_numpy()
        j = adata.var_names.get_loc(gene)
        for group, row in result[result.names == gene].set_index("group").iterrows():
            inside = cells & (adata.obs.cluster == group).to_numpy()
            ref = mannwhitneyu(x[inside, j], x[cells & ~inside, j], method="asymptotic", use_continuity=False)
            assert row.U == ref.statistic
            assert row.pvals == pytest.approx(ref.pvalue, rel=1e-9)
            assert (row.n_measured_group, row.n_measured_rest) == (inside.sum(), (cells & ~inside).sum())
    # one BH family per group over both routes
    assert result.groupby("group").size().eq(6).all()
    assert (result.groupby("group").scores.apply(lambda s: s.is_monotonic_decreasing)).all()


def test_groups_with_too_few_measuring_cells_are_unassessed():
    adata, fm = _adata()
    in_a = (adata.obs.cluster == "a").to_numpy()
    adata.obs.loc[in_a, "sample"] = "s2"
    adata.obs.loc[adata.obs.index[in_a][:2], "sample"] = "s3"  # 2 < min_measured cells of "a" measure g5
    result = run_markers_scanpy(adata, groupby="cluster", layer="X", measurement=fm)
    assert set(result.loc[result.names == "g5", "group"]) == {"b", "c"}
    assert len(result[result.names == "g4"]) == 3


def test_log10_pval_is_finite_where_pval_underflows():
    n = 4000
    x = np.arange(1, n + 1) / n * 5.0  # "lo" cells all below "hi" cells
    adata = AnnData(np.c_[x, x[::-1]], obs=pd.DataFrame({"cluster": pd.Categorical(["lo"] * (n // 2) + ["hi"] * (n // 2))},
                                                      index=[f"c{i}" for i in range(n)]),
                    var=pd.DataFrame(index=["up", "down"]))
    result = run_markers_scanpy(adata, groupby="cluster", layer="X")
    assert (result.pvals == 0).all()
    assert np.isfinite(result.log10_pval).all() and (result.log10_pval < -300).all()
    assert np.isfinite(result.log10_pval_adj).all()
    assert set(result.U) == {0.0, (n // 2) ** 2}


def test_evidence_command_discovers_registered_measurement(tmp_path, monkeypatch):
    from scoop.commands import annotation_evidence

    run = tmp_path / "runs/test"
    (run / "config").mkdir(parents=True)
    config = run / "config/run.yaml"
    config.write_text(yaml.safe_dump({"run": {"run_id": "test"}, "annotation": {"decision": {"cluster_key": "cluster"}}}))
    adata, fm = _adata()
    adata.layers["log1p_norm"] = adata.X.copy()
    plain_input = run / "artifacts/plain.h5ad"
    plain_input.parent.mkdir()
    adata.write_h5ad(plain_input)
    monkeypatch.setattr(sys, "argv", ["evidence", "--config", str(config), "--input", str(plain_input)])
    assert annotation_evidence.main()["measurement"] == {"status": "not_declared"}

    path = fm.save(run / "artifacts/feature_measurement.npz")
    register_file(adata, key="feature_measurement", path=path, category="artifacts", schema=FORMAT)
    adata.write_h5ad(run / "artifacts/adata_programs.h5ad")
    monkeypatch.setattr(sys, "argv", ["evidence", "--config", str(config)])
    result = annotation_evidence.main()
    assert result["measurement"]["status"] == "applied" and result["measurement"]["n_genes_masked"] == 2
    assert str(path.resolve()) in result["inputs"]
    markers = pd.read_parquet(run / "annotation/tables/cluster_markers.parquet")
    g5 = markers[markers.names == "g5"]
    assert (g5.n_measured_group + g5.n_measured_rest == (adata.obs["sample"] == "s3").sum()).all()
    assert "feature_measurement" in read_h5ad(run / "artifacts/adata_annotation_evidence.h5ad").uns["file_registry"]["artifacts"]
