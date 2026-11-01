from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from scsp_agent_sop.tumor_cnv import prepare_tumor_fastcnv_reference, run_tumor_fastcnv_workflow


def test_prepare_tumor_fastcnv_reference_marks_normal_pools() -> None:
    obs = pd.DataFrame(
        {
            "sample_id": ["s1"] * 4 + ["s2"] * 4,
            "cell_type_lvl1": ["T", "Myeloid", "Epithelial", "Epithelial", "B", "Endothelial", "Epithelial", "Epithelial"],
            "annotation_confidence": ["high"] * 8,
            "tumor_status": ["normal", "normal", "normal", "tumor", "normal", "normal", "mixed", "tumor"],
        },
        index=[f"c{i}" for i in range(8)],
    )
    adata = ad.AnnData(X=np.ones((8, 3)), obs=obs, var=pd.DataFrame(index=["G1", "G2", "G3"]))

    plan = prepare_tumor_fastcnv_reference(
        adata,
        major_lineage_key="cell_type_lvl1",
        sample_key="sample_id",
        parenchymal_lineages=("Epithelial",),
        normal_status_key="tumor_status",
        min_reference_cells=2,
    )

    assert plan.should_run_fastcnv
    assert plan.reference_labels == ["normal_nonparenchymal", "normal_parenchymal"]
    assert plan.n_normal_nonparenchymal == 4
    assert plan.n_normal_parenchymal == 1
    assert plan.n_candidate_cells == 3
    assert adata.obs.loc["c0", "fastcnv_reference_pool"] == "normal_nonparenchymal"
    assert adata.obs.loc["c2", "fastcnv_reference_pool"] == "normal_parenchymal"
    assert adata.obs.loc["c3", "fastcnv_reference_pool"] == "candidate_tumor_or_mixed"


def test_unresolved_low_confidence_and_filtered_cells_are_never_reference():
    import pytest
    obj = ad.AnnData(np.ones((5, 3)), obs=pd.DataFrame({
        'sample_id': ['s'] * 5,
        'cell_type_lvl1': ['Unresolved', 'T', 'T', 'T', 'Epithelial'],
        'annotation_confidence': ['high', 'low', 'high', 'high', 'high'],
        'final_use': [True, True, False, True, True],
    }, index=[str(i) for i in range(5)]))
    plan = prepare_tumor_fastcnv_reference(obj, min_reference_cells=1)
    assert plan.n_reference_cells == 1
    assert obj.obs.fastcnv_reference_pool.eq('normal_nonparenchymal').tolist() == [False, False, False, True, False]
    del obj.obs['annotation_confidence']
    with pytest.raises(ValueError, match='reviewed'):
        prepare_tumor_fastcnv_reference(obj)


def test_committed_and_legacy_statuses_respect_explicit_reference_review():
    obs = pd.DataFrame({
        'sample_id': ['s'] * 5,
        'cell_type_lvl1': ['T'] * 5,
        'annotation_status': ['accepted', 'accepted', 'labeled', 'review_required', 'unlabeled'],
        'reviewed_reference': [True, False, True, True, True],
    }, index=[f'c{i}' for i in range(5)])
    obj = ad.AnnData(np.ones((5, 3)), obs=obs)
    plan = prepare_tumor_fastcnv_reference(
        obj, reference_eligible_key='reviewed_reference', min_reference_cells=2,
    )
    assert plan.n_reference_cells == 2
    assert obj.obs.fastcnv_reference_pool.eq('normal_nonparenchymal').tolist() == [True, False, True, False, False]


def test_duplicate_symbols_cannot_silently_select_a_genomic_column(tmp_path):
    import pytest
    obj = ad.AnnData(np.array([[1, 100, 2], [3, 200, 4]]),
        var=pd.DataFrame(index=['ATRIP', 'ATRIP', 'G2']))
    with pytest.raises(ValueError, match='CNV genomic positions require unique'):
        run_tumor_fastcnv_workflow(obj, pd.DataFrame(), output_dir=tmp_path)
