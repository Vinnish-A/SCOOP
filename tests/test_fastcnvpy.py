from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

from fastcnvpy import FastCNVConfig, run_fastcnv, run_fastcnv_anndata, run_fastcnv_pooled_anndata
from fastcnvpy.core import _compute_genomic_scores


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
