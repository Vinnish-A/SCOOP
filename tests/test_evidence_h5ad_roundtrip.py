import json

import anndata as ad
import numpy as np
import pandas as pd

from scsp_agent_sop.annotation_decision.evidence_bundle import build_evidence_bundle


def test_h5ad_evidence_remains_json_serializable(tmp_path):
    a = ad.AnnData(np.ones((3, 2)), obs=pd.DataFrame(
        {"cluster_identity": ["0"] * 3}, index=["a", "b", "c"]))
    a.uns["file_registry"] = {"tables": {"markers": {"bytes": 123, "columns": ["gene", "score"]}}}
    a.uns["top_markers_0"] = ["CLEC9A", "XCR1"]
    a.write_h5ad(tmp_path / "input.h5ad")
    bundle = build_evidence_bundle(ad.read_h5ad(tmp_path / "input.h5ad"), run_id="r",
        organism="human", tissue="PBMC", disease_context=None, is_tumor=False,
        cluster_key="cluster_identity")
    payload = json.loads(bundle.to_json())
    assert payload["clusters"][0]["top_markers"] == ["CLEC9A", "XCR1"]
    assert payload["file_registry"]["tables"]["markers"]["bytes"] == 123
