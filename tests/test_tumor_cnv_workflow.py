from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from scsp_agent_sop.tumor_cnv import prepare_tumor_fastcnv_reference


def test_prepare_tumor_fastcnv_reference_marks_normal_pools() -> None:
    obs = pd.DataFrame(
        {
            "sample_id": ["s1"] * 4 + ["s2"] * 4,
            "cell_type_lvl1": ["T", "Myeloid", "Epithelial", "Epithelial", "B", "Endothelial", "Epithelial", "Epithelial"],
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
