import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scsp_agent_sop.annotation import run_markers_scanpy
from scsp_agent_sop.annotation_decision.evidence_bundle import build_evidence_bundle


def test_differential_export_uses_all_requested_layer_genes_and_ignores_raw_axis():
    matrix = np.log1p(np.array([[6, 5, 4, 0, 0, 0]] * 4 + [[0, 0, 0, 4, 5, 6]] * 4, dtype=float))
    obj = AnnData(matrix, obs=pd.DataFrame({"cluster_identity": pd.Categorical(["0"] * 4 + ["1"] * 4)},
                    index=[f"c{i}" for i in range(8)]), var=pd.DataFrame(index=[f"g{i}" for i in range(6)]))
    raw = AnnData(np.ones((8, 1)), obs=obj.obs.copy(), var=pd.DataFrame(index=["unrelated_raw_gene"]))
    obj.raw = raw
    obj.layers["log1p_norm"] = matrix.copy()
    before = obj.X.copy()
    result = run_markers_scanpy(obj)
    assert len(result) == 12 and set(result.names) == set(obj.var_names)
    assert result.scores.notna().all() and result.pvals_adj.notna().all()
    np.testing.assert_array_equal(obj.X, before)
    with pytest.raises(ValueError, match="layer is missing"):
        run_markers_scanpy(obj, layer="missing")


def test_registered_cluster_deg_reaches_review_template_bundle():
    obj = AnnData(np.ones((4, 2)), obs=pd.DataFrame({"cluster_identity": ["0"] * 4}, index=list("abcd")))
    obj.uns["file_registry"] = {"tables": {"cluster_markers": {"schema": "cluster_markers.v1", "path": "genes.parquet"},
                                          "scores": {"schema": "marker_means.v1", "path": "scores.parquet"}}}
    bundle = build_evidence_bundle(obj, run_id="r", organism="human", tissue="brain", disease_context=None,
                                   is_tumor=False, cluster_key="cluster_identity")
    assert bundle.clusters[0].marker_refs == ("cluster_markers",)
