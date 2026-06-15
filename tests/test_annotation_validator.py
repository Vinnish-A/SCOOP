from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from scsp_agent_sop.annotation_decision.committer import commit_annotation_decisions
from scsp_agent_sop.annotation_decision.decision_schema import AnnotationDecision
from scsp_agent_sop.annotation_decision.evidence_bundle import build_evidence_bundle
from scsp_agent_sop.annotation_decision.validator import validate_annotation_decision


def _decision(**overrides) -> AnnotationDecision:
    data = {
        "schema_version": "scoop.annotation_decision.v1",
        "run_id": "run1",
        "cluster_id": "0",
        "cluster_key": "cluster_identity",
        "parent_label": "Epithelial",
        "canonical_label": "Malignant epithelial cell",
        "cell_state": None,
        "functional_modifier": None,
        "final_label": "Malignant epithelial cell",
        "confidence": "medium",
        "evidence_refs": {},
        "positive_markers": ("EPCAM",),
        "negative_markers_absent": (),
        "conflicts": (),
        "review_required": False,
        "reason": "synthetic",
    }
    data.update(overrides)
    return AnnotationDecision.from_dict(data)


def _tumor_bundle(has_cnv: bool):
    adata = AnnData(np.ones((4, 2)))
    obs = {"cluster_identity": ["0", "0", "1", "1"], "sample_id": ["s1", "s2", "s1", "s2"]}
    if has_cnv:
        obs["fastcnv_tumor_evidence"] = ["tumor", "tumor", "normal", "normal"]
        obs["fastcnv_cnv_fraction"] = [0.5, 0.4, 0.01, 0.02]
    adata.obs = pd.DataFrame(obs, index=[f"c{i}" for i in range(4)])
    return adata, build_evidence_bundle(
        adata,
        run_id="run1",
        organism="human",
        tissue="brain",
        disease_context="glioblastoma",
        is_tumor=True,
        cluster_key="cluster_identity",
    )


def test_tumor_malignant_label_without_cnv_evidence_forces_review() -> None:
    _, bundle = _tumor_bundle(False)
    result = validate_annotation_decision(_decision(), bundle)
    assert result.valid is True
    assert result.review_required is True
    assert result.accepted is False


def test_tumor_malignant_label_with_cnv_evidence_can_pass() -> None:
    _, bundle = _tumor_bundle(True)
    result = validate_annotation_decision(_decision(evidence_refs={"cnv": ("fastcnv_tumor_evidence",)}), bundle)
    assert result.valid is True
    assert result.review_required is False
    assert result.accepted is True


def test_marker_conflict_prevents_high_confidence() -> None:
    adata, bundle = _tumor_bundle(True)
    del adata
    decision = _decision(confidence="high", conflicts=("PTPRC conflict",), evidence_refs={"cnv": ("fastcnv_tumor_evidence",)})
    result = validate_annotation_decision(decision, bundle)
    assert result.review_required is True
    assert result.accepted is False


def test_committer_writes_expected_obs_fields(tmp_path) -> None:
    adata, bundle = _tumor_bundle(True)
    decision = _decision(evidence_refs={"cnv": ("fastcnv_tumor_evidence",)})
    commit_annotation_decisions(adata, [decision], bundle, run_root=tmp_path, output_dir=tmp_path / "tables")
    for field in ["cell_type_lvl1", "cell_type_lvl2", "cell_type_lvl3", "cell_state", "annotation_confidence", "annotation_status"]:
        assert field in adata.obs
    mask = adata.obs["cluster_identity"].astype(str) == "0"
    assert set(adata.obs.loc[mask, "cell_type_lvl3"]) == {"Malignant epithelial cell"}
    assert set(adata.obs.loc[mask, "annotation_status"]) == {"accepted"}
    assert (tmp_path / "tables" / "annotation_decisions.tsv").exists()
    assert (tmp_path / "logs" / "decision_log.jsonl").exists()
