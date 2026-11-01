from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
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
    data["evidence_refs"] = {"deg": ("manual:synthetic_cluster_DEG_fixture",), **data["evidence_refs"]}
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


# Phase G C: lineage markers decide malignancy and CNV is supporting only, so
# the former "without CNV forces review" test is replaced by these two.
def test_tumor_label_without_declared_malignancy_forces_review() -> None:
    _, bundle = _tumor_bundle(False)
    result = validate_annotation_decision(_decision(), bundle)
    assert result.valid is True
    assert result.review_required is True
    assert result.accepted is False
    assert any("malignancy_status" in reason for reason in result.reasons)


def test_declared_lineage_malignancy_passes_without_cnv() -> None:
    _, bundle = _tumor_bundle(False)
    result = validate_annotation_decision(
        _decision(malignancy_status="malignant", malignancy_basis=("lineage_markers",)), bundle)
    assert result.accepted is True


def test_tumor_malignant_label_with_cnv_evidence_can_pass() -> None:
    _, bundle = _tumor_bundle(True)
    result = validate_annotation_decision(_decision(evidence_refs={"cnv": ("fastcnv_tumor_evidence",)},
                                                    malignancy_status="malignant", malignancy_basis=("lineage_markers",)), bundle)
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
    decision = _decision(evidence_refs={"cnv": ("fastcnv_tumor_evidence",)},
                         malignancy_status="malignant", malignancy_basis=("lineage_markers",))
    commit_annotation_decisions(adata, [decision], bundle, run_root=tmp_path, output_dir=tmp_path / "tables")
    for field in ["cell_type_lvl1", "cell_type_lvl2", "cell_type_lvl3", "cell_state", "annotation_confidence", "annotation_status"]:
        assert field in adata.obs
    mask = adata.obs["cluster_identity"].astype(str) == "0"
    assert set(adata.obs.loc[mask, "cell_type_lvl3"]) == {"Malignant epithelial cell"}
    assert set(adata.obs.loc[mask, "annotation_status"]) == {"accepted"}
    assert (tmp_path / "tables" / "annotation_decisions.tsv").exists()
    assert (tmp_path / "logs" / "decision_log.jsonl").exists()
    assert adata.uns["annotation_review"]["reviewed"] == 2
    assert adata.uns["annotation_review"]["unreviewed"] == 2


def test_rejected_replacement_revokes_old_accepted_reference_status(tmp_path):
    adata, bundle = _tumor_bundle(False)
    adata.obs["annotation_status"] = pd.Categorical(["accepted"] * 4)
    adata.obs["cell_type_lvl3"] = pd.Categorical(["old"] * 4)
    commit_annotation_decisions(adata, [_decision()], bundle, run_root=tmp_path, output_dir=tmp_path / "tables")
    assert adata.obs.annotation_status.tolist() == ["review_required", "review_required", "accepted", "accepted"]
    assert adata.obs.annotation_reviewed.tolist() == [True, True, False, False]
    assert adata.uns["annotation_review"]["accepted"] == 0


def test_overlapping_decisions_fail_before_obs_mutation(tmp_path):
    adata, bundle = _tumor_bundle(True)
    before = adata.obs.copy(deep=True)
    with pytest.raises(ValueError, match="overlap"):
        commit_annotation_decisions(adata, [_decision(), _decision()], bundle, run_root=tmp_path, output_dir=tmp_path / "tables")
    pd.testing.assert_frame_equal(before, adata.obs)


def test_broad_identity_does_not_accept_fine_identity_or_assign_normality(tmp_path):
    adata, bundle = _tumor_bundle(False)
    decision = _decision(parent_label="Glial_RNA", canonical_label="Glial_RNA", final_label="Glial_Unresolved", identity_scope="broad")
    commit_annotation_decisions(adata, [decision], bundle, run_root=tmp_path, output_dir=tmp_path / "tables")
    assert adata.obs.broad_annotation_status.tolist()[:2] == ["accepted"] * 2
    assert adata.obs.fine_annotation_status.tolist()[:2] == ["review_required"] * 2
    assert "malignancy_status" not in adata.obs
    assert adata.uns["annotation_review"]["fine_accepted"] == 0


@pytest.mark.parametrize("scope,expected_fine", [(None, "review_required"), ("fine", "accepted")])
def test_fine_acceptance_requires_explicit_scope(tmp_path, scope, expected_fine):
    adata, bundle = _tumor_bundle(False)
    fields = dict(parent_label="T", canonical_label="T", final_label="T_cell")
    if scope is not None:
        fields["identity_scope"] = scope
    decision = _decision(**fields)
    commit_annotation_decisions(adata, [decision], bundle, run_root=tmp_path, output_dir=tmp_path / "tables")
    assert adata.obs.annotation_status.tolist()[:2] == ["accepted"] * 2
    assert adata.obs.broad_annotation_status.tolist()[:2] == ["accepted"] * 2
    assert adata.obs.fine_annotation_status.tolist()[:2] == [expected_fine] * 2


def test_marker_scores_without_differential_evidence_require_review():
    _, bundle = _tumor_bundle(False)
    data = _decision(parent_label="Mast", canonical_label="Mast", final_label="Mast").to_dict()
    data["evidence_refs"] = {"scores": ["manual:KIT_mean_score"]}
    result = validate_annotation_decision(AnnotationDecision.from_dict(data), bundle)
    assert not result.accepted
    assert any("differential-gene" in reason for reason in result.reasons)


def test_explicit_unresolved_decision_records_candidate_without_accepting_it(tmp_path):
    adata, bundle = _tumor_bundle(False)
    decision = _decision(parent_label="Lymphoid", canonical_label="B_RNA", final_label="B_lineage_Unresolved",
                         identity_scope="unresolved", confidence="low", review_required=True)
    commit_annotation_decisions(adata, [decision], bundle, run_root=tmp_path, output_dir=tmp_path / "tables")
    assert adata.obs.cell_type_lvl3.tolist()[:2] == ["B_lineage_Unresolved"] * 2
    assert adata.obs.annotation_status.tolist()[:2] == ["review_required"] * 2
    assert adata.uns["annotation_review"]["accepted"] == 0
