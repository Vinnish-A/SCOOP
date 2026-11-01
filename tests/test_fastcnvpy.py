from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad
import pytest
import scipy.sparse as sp

from fastcnvpy import FastCNVConfig, run_fastcnv, run_fastcnv_anndata, run_fastcnv_pooled_anndata
from fastcnvpy.core import (
    _compute_genomic_scores,
    _select_genes_like_fastcnv,
    classify_arm_cnv,
    cnv_per_chromosome_arm,
    merge_cnv_clusters,
)
from fastcnvpy.reference import normalize_gene_metadata


def _inputs():
    rng = np.random.default_rng(3)
    genes = [f"G{i:04d}" for i in range(360)]
    cells = [f"cell{i:03d}" for i in range(30)]
    counts = rng.poisson(5, size=(len(genes), len(cells))).astype(float)
    counts[:80, 15:] *= 2.0
    counts = pd.DataFrame(np.round(counts), index=genes, columns=cells)
    chroms = np.repeat([1, 2, 3, 4], 90)
    arms = np.tile(np.repeat(["p", "q"], 45), 4)
    starts = np.tile(np.arange(90) * 1000 + 1, 4)
    gene_metadata = pd.DataFrame(
        {
            "hgnc_symbol": genes,
            "chromosome_name": chroms.astype(str),
            "start_position": starts,
            "end_position": starts + 100,
            "gene_biotype": "protein_coding",
            "chr_arm": arms,
        }
    )
    obs = pd.DataFrame({"annot": ["normal"] * 15 + ["tumor"] * 15}, index=cells)
    return counts, gene_metadata, obs


def test_fastcnvpy_outputs_fastcnv_core_slots() -> None:
    counts, gene_metadata, obs = _inputs()

    result = run_fastcnv(
        counts,
        gene_metadata,
        obs=obs,
        reference_var="annot",
        reference_label="normal",
        config=FastCNVConfig(window_size=20, window_step=5, top_n_genes=120, cluster_k=2),
        sample_name="tiny",
    )

    assert result.raw_genomic_scores.shape == result.genomic_scores.shape
    assert result.raw_genomic_scores.shape[1] == counts.shape[1]
    assert result.cell_metadata.index.tolist() == counts.columns.tolist()
    assert "cnv_fraction" in result.cell_metadata.columns
    assert "1.p_CNV" in result.cell_metadata.columns
    assert "X.q_CNV" in result.cell_metadata.columns
    assert "cnv_clusters" in result.cell_metadata.columns
    assert result.arm_cnv.shape == (46, counts.shape[1])
    assert set(result.window_metadata.columns) == {"window", "chrom_arm", "n_genes", "genes"}
    assert result.manifest["n_windows"] == result.raw_genomic_scores.shape[0]


def test_genomic_scores_match_r_duplicate_rowname_lookup() -> None:
    norm_counts = pd.DataFrame(
        [[1.0, 3.0], [10.0, 30.0], [5.0, 7.0]],
        index=["G1", "G1", "G2"],
        columns=["cell1", "cell2"],
    )

    scores = _compute_genomic_scores(norm_counts, {"1.p1": ["G1", "G2"], "1.p2": ["G1"]})

    np.testing.assert_allclose(scores[:, 0], [3.0, 5.0])
    np.testing.assert_allclose(scores[:, 1], [1.0, 3.0])


def test_anndata_path_matches_dataframe_path() -> None:
    counts, gene_metadata, obs = _inputs()
    cfg = FastCNVConfig(window_size=20, window_step=5, top_n_genes=120, cluster_k=2)
    frame_result = run_fastcnv(
        counts,
        gene_metadata,
        obs=obs,
        reference_var="annot",
        reference_label="normal",
        config=cfg,
        sample_name="frame",
        compute_clusters=False,
    )
    adata = ad.AnnData(
        X=sp.csc_matrix(counts.T.to_numpy()),
        obs=obs.copy(),
        var=pd.DataFrame(index=counts.index),
    )

    adata_result = run_fastcnv_anndata(
        adata,
        gene_metadata,
        reference_var="annot",
        reference_label="normal",
        config=cfg,
        sample_name="adata",
        compute_clusters=False,
    )

    pd.testing.assert_frame_equal(adata_result.raw_genomic_scores, frame_result.raw_genomic_scores)
    pd.testing.assert_frame_equal(adata_result.genomic_scores, frame_result.genomic_scores)
    pd.testing.assert_frame_equal(adata_result.arm_cnv, frame_result.arm_cnv)


def test_pooled_anndata_splits_samples_with_pooled_reference() -> None:
    counts, gene_metadata, obs = _inputs()
    obs = obs.copy()
    obs["sample_id"] = ["s1"] * 15 + ["s2"] * 15
    adata = ad.AnnData(
        X=sp.csc_matrix(counts.T.to_numpy()),
        obs=obs,
        var=pd.DataFrame(index=counts.index),
    )

    result = run_fastcnv_pooled_anndata(
        adata,
        gene_metadata,
        sample_key="sample_id",
        reference_var="annot",
        reference_label="normal",
        config=FastCNVConfig(window_size=20, window_step=5, top_n_genes=120, cluster_k=2),
        compute_clusters=False,
    )

    assert result.manifest["schema_version"] == "fastcnvpy.pooled_result.v1"
    assert result.manifest["samples"] == ["s1", "s2"]
    assert set(result.per_sample) == {"s1", "s2"}
    assert result.cell_metadata.index.tolist() == adata.obs_names.tolist()
    assert "cnv_fraction" in result.cell_metadata.columns
    assert result.arm_cnv.shape[0] == 46
    assert result.window_metadata.shape[0] > 0


def test_qualified_arms_and_blocked_pooled_scores_match():
    counts, metadata, obs = _inputs()
    obs['sample_id'] = ['ref'] * 15 + ['candidate'] * 15
    obj = ad.AnnData(sp.csr_matrix(counts.T.to_numpy()), obs=obs, var=pd.DataFrame(index=counts.index))
    cfg = FastCNVConfig(window_size=20, window_step=5, top_n_genes=120)
    kwargs = dict(reference_var='annot', reference_label='normal', config=cfg, compute_clusters=False)
    baseline = run_fastcnv_pooled_anndata(obj, metadata, cell_block_size=100, **kwargs)
    qualified = metadata.copy()
    qualified['chr_arm'] = qualified.chromosome_name + '.' + qualified.chr_arm
    blocked = run_fastcnv_pooled_anndata(obj, qualified, cell_block_size=3, **kwargs)
    np.testing.assert_allclose(baseline.arm_cnv, blocked.arm_cnv)
    obj.X[15:, :80] *= 5
    changed = run_fastcnv_pooled_anndata(obj, metadata, cell_block_size=3, **kwargs)
    np.testing.assert_allclose(baseline.per_sample['ref'].raw_genomic_scores, changed.per_sample['ref'].raw_genomic_scores)
    assert 'cnv_clusters' not in blocked.cell_metadata


def test_windowless_arm_is_unassessed_and_does_not_block_merging() -> None:
    scores = pd.DataFrame([[0.0, 0.2], [-0.3, 0.0], [0.4, 0.4]], index=["1.p1", "1.q1", "13.q1"], columns=["c1", "c2"])

    arm = cnv_per_chromosome_arm(scores)
    classes = classify_arm_cnv(arm, (-0.1, 0.0, 0.1))

    assert arm.loc["13.p_CNV"].isna().all()
    assert classes.loc["13.p_CNV_classification"].tolist() == ["unassessed", "unassessed"]
    assert classes.loc["13.q_CNV_classification"].tolist() == ["gain", "gain"]
    meta = pd.DataFrame(
        {"cnv_clusters": [1, 1, 2, 2], "1.p_CNV": [0.2, 0.21, 0.19, 0.2], "1.q_CNV": [-0.3, -0.31, -0.29, -0.3], "13.p_CNV": np.nan},
        index=list("abcd"),
    )
    assert merge_cnv_clusters(meta).tolist() == [1, 1, 1, 1]


def test_gene_selection_keeps_fastcnv_gene_set_without_repeated_positions() -> None:
    rows = [(f"g{c}{a}{i}", c, i * 1000 + (0 if a == "p" else 10**8), i * 1000 + 500 + (0 if a == "p" else 10**8), a)
            for c in (1, 2) for a in "pq" for i in range(150)]
    metadata = pd.DataFrame(rows, columns=["hgnc_symbol", "chromosome_num", "start_position", "end_position", "chr_arm"])
    genes = metadata["hgnc_symbol"].tolist()
    average = pd.Series(np.random.default_rng(0).gamma(2.0, 2.0, len(genes)), index=genes)

    selected = _select_genes_like_fastcnv(common_genes=genes, average_expression=average, gene_metadata=metadata,
                                          top_n_genes=100, genes_to_force=None, chr_arms_to_force=None, region_to_force=None)

    final = selected.final_selected_genes
    assert len(final) == len(set(final))
    # fastCNV's set: the top-N genes plus the genome-wide top 200 added by the bare p/q loop.
    expected = set(average.sort_values(ascending=False).index[:200])
    assert set(final) == expected | set(selected.top_expr_genes)


def test_normalized_metadata_applies_biomart_scope_and_rejects_repeated_rows() -> None:
    frame = pd.DataFrame({"hgnc_symbol": ["A", "XG", "YGENE", "MT-X", ""], "chromosome_num": [1, 23, 24, 25, 2],
                          "start_position": [10, 10, 5, 1, 3], "end_position": [20, 20, 9, 2, 4], "chr_arm": ["p", "p", "p", "p", "q"]})

    assert normalize_gene_metadata(frame)["hgnc_symbol"].tolist() == ["A", "XG"]
    with pytest.raises(ValueError, match="repeats 1"):
        normalize_gene_metadata(pd.concat([frame, frame.iloc[[0]].assign(chr_arm="1.p")]))


def test_anndata_refuses_log_normalised_x_as_counts() -> None:
    counts, gene_metadata, obs = _inputs()
    adata = ad.AnnData(X=sp.csr_matrix(np.log1p(counts.T.to_numpy())), obs=obs.copy(), var=pd.DataFrame(index=counts.index))
    kwargs = dict(reference_var="annot", reference_label="normal", config=FastCNVConfig(window_size=20, window_step=5, top_n_genes=120),
                  compute_clusters=False)

    with pytest.raises(ValueError, match="no 'counts' layer"):
        run_fastcnv_anndata(adata, gene_metadata, **kwargs)
    adata.layers["counts"] = sp.csr_matrix(counts.T.to_numpy())
    assert run_fastcnv_anndata(adata, gene_metadata, **kwargs).raw_genomic_scores.shape[1] == adata.n_obs
