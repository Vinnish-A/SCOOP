"""Phase G annotation contracts: content-bound evidence, declared malignancy, outcome schema."""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd
import pytest
import yaml
from anndata import AnnData, read_h5ad

from scoop.commands import annotation_commit, annotation_review
from scoop.commands.deliver import review_coverage
from scsp_agent_sop.annotation_decision.committer import commit_annotation_decisions
from scsp_agent_sop.annotation_decision.decision_schema import AnnotationDecision
from scsp_agent_sop.annotation_decision.evidence_bundle import AnnotationEvidenceBundle, build_evidence_bundle
from scsp_agent_sop.annotation_decision.validator import validate_annotation_decision
from scsp_agent_sop.storage import register_file


def _adata(tmp_path) -> AnnData:
    obs = pd.DataFrame({"cluster_identity": ["0"] * 3 + ["1"] * 3,
                        "sample_key": ["s1", "s1", "s2", "s1", "s2", "s3"]}, index=[f"c{i}" for i in range(6)])
    adata = AnnData(np.ones((6, 2)), obs=obs)
    pd.DataFrame({"group": ["0", "0", "0"], "names": ["CD2", "CD3E", "LYZ"], "scores": [5.0, 9.0, 1.0]}).to_parquet(
        tmp_path / "deg.parquet")
    register_file(adata, key="cluster_markers", path=tmp_path / "deg.parquet", schema="cluster_markers.v1",
                  extra={"label_key": "cluster_identity"})
    # DEG of another partition that happens to share group names with this one.
    pd.DataFrame({"group": ["1", "7"], "names": ["MS4A1", "PDGFRB"], "scores": [3.0, 2.0]}).to_csv(
        tmp_path / "parent.tsv", sep="\t", index=False)
    register_file(adata, key="parent_markers", path=tmp_path / "parent.tsv", schema="subcluster_markers.v1",
                  extra={"label_key": "parent_cluster"})
    pd.DataFrame({"obs_name": ["c3"], "doublet_score": [0.9]}).to_csv(tmp_path / "qc.tsv", sep="\t", index=False)
    register_file(adata, key="doublet_scores", path=tmp_path / "qc.tsv", schema="doublet_scores.v1")
    (tmp_path / "cnv.json").write_text("{}")
    register_file(adata, key="cnv_summary", path=tmp_path / "cnv.json", schema="cnv_summary.supporting_only.v1")
    return adata


def _bundle(adata, **kwargs) -> AnnotationEvidenceBundle:
    return build_evidence_bundle(adata, run_id="run1", organism="human", tissue="brain", disease_context=None,
                                 is_tumor=False, cluster_key="cluster_identity", sample_key="sample_key", **kwargs)


def _decision(**overrides) -> AnnotationDecision:
    data = {"schema_version": "scoop.annotation_decision.v1", "run_id": "run1", "cluster_id": "0",
            "cluster_key": "cluster_identity", "parent_label": "Lymphoid", "canonical_label": "T", "cell_state": None,
            "functional_modifier": None, "final_label": "T_cell", "confidence": "medium",
            "evidence_refs": {"deg": ["cluster_markers"]}, "positive_markers": ["CD3E"], "negative_markers_absent": [],
            "conflicts": [], "review_required": False, "reason": "synthetic"}
    data.update(overrides)
    return AnnotationDecision.from_dict(data)


def test_deg_attaches_only_to_units_its_table_tests(tmp_path):
    bundle = _bundle(_adata(tmp_path))
    clusters = bundle.cluster_by_id()
    assert bundle.deg_tables["cluster_markers"] == {"label_key": "cluster_identity", "groups": ["0"],
                                                    "sha256": bundle.deg_tables["cluster_markers"]["sha256"]}
    assert clusters["0"].deg_refs == ("cluster_markers",) and clusters["0"].marker_refs == ("cluster_markers",)
    # Unit 1 shares a group name with parent_markers but that table tests another partition.
    assert clusters["1"].deg_refs == () and clusters["1"].marker_refs == ()
    assert clusters["0"].top_markers == ("CD3E", "CD2", "LYZ")
    assert bundle.partition_sha256 and bundle.schema_version.endswith("v2")


def test_deg_citation_must_resolve_to_a_tested_group(tmp_path):
    bundle = _bundle(_adata(tmp_path))
    uncovered = validate_annotation_decision(_decision(cluster_id="1"), bundle)
    assert not uncovered.accepted and any("does not test this unit" in r for r in uncovered.reasons)
    other_partition = validate_annotation_decision(_decision(cluster_id="1", evidence_refs={"deg": ["parent_markers"]}), bundle)
    assert not other_partition.accepted
    parent_group = validate_annotation_decision(
        _decision(cluster_id="1", evidence_refs={"deg": ["parent_markers::7"]}, parent_unit="parent_cluster:7"), bundle)
    assert parent_group.accepted, parent_group.reasons
    assert not validate_annotation_decision(_decision(cluster_id="1", evidence_refs={"deg": ["parent_markers::9"]}), bundle).accepted


def test_sample_support_uses_configured_key_and_per_sample_minimum(tmp_path):
    adata = _adata(tmp_path)
    cluster = _bundle(adata, min_cells_per_sample=2).cluster_by_id()["0"]
    assert (cluster.n_samples, cluster.n_samples_raw) == (1, 2)
    unknown = build_evidence_bundle(adata, run_id="run1", organism="human", tissue="brain", disease_context=None,
                                    is_tumor=False, cluster_key="cluster_identity")
    result = validate_annotation_decision(_decision(confidence="high"), unknown, min_cells_high_confidence=1)
    assert any("sample support is unknown" in r for r in result.reasons)


def test_commit_refuses_stale_partition_and_count_mismatch(tmp_path):
    adata = _adata(tmp_path)
    bundle = _bundle(adata)
    moved = adata.copy()
    moved.obs.loc["c2", "cluster_identity"] = "1"
    with pytest.raises(ValueError, match="different partition"):
        commit_annotation_decisions(moved, [_decision()], bundle, run_root=tmp_path, output_dir=tmp_path / "t")
    legacy = bundle.to_dict()
    for key in ("partition_sha256", "deg_tables", "composition_keys"):
        legacy.pop(key)
    legacy = AnnotationEvidenceBundle.from_dict(legacy)
    assert legacy.partition_sha256 is None
    with pytest.raises(ValueError, match="has 2 cells but the evidence bundle has 3"):
        commit_annotation_decisions(moved, [_decision()], legacy, run_root=tmp_path, output_dir=tmp_path / "t")
    assert "annotation_status" not in moved.obs


def test_malignancy_is_declared_and_lineage_based(tmp_path):
    bundle = _bundle(_adata(tmp_path))
    malignant = dict(parent_label="Malignant", canonical_label="GBM_Malignant", final_label="GBM_Malignant_AC")

    def check(**fields):
        return validate_annotation_decision(_decision(**{**malignant, **fields}), bundle)

    assert check(malignancy_status="malignant", malignancy_basis=["lineage_markers"]).accepted
    assert not check().accepted
    assert not check(malignancy_status="normal").accepted
    assert not check(malignancy_status="malignant", malignancy_basis=["cnv"]).accepted
    assert not check(malignancy_status="malignant", malignancy_basis=["lineage_markers"], evidence_refs={}).accepted
    assert not check(malignancy_status="malignant", malignancy_basis=["lineage_markers", "cnv"],
                     evidence_refs={"deg": ["cluster_markers"], "cnv": ["fastcnv_tumor_evidence"]}).accepted
    assert check(malignancy_status="malignant", malignancy_basis=["lineage_markers", "cnv"],
                 evidence_refs={"deg": ["cluster_markers"], "cnv": ["cnv_summary"]}).accepted
    assert not check(malignancy_status="malignant-ish", malignancy_basis=["lineage_markers"]).accepted
    # A declared status gates labels without tumour words (C1-10).
    assert not validate_annotation_decision(_decision(malignancy_status="probable_malignant"), bundle).accepted
    assert validate_annotation_decision(_decision(**malignant), bundle, tumor_gate={"enabled": False}).accepted


def test_technical_exclusion_is_a_reviewed_outcome(tmp_path):
    adata = _adata(tmp_path)
    bundle = _bundle(adata)
    technical = _decision(cluster_id="1", parent_label="Technical", canonical_label="Technical", final_label="Doublet_T_Myeloid",
                          disposition="technical_exclusion", technical_kind="doublet",
                          evidence_refs={"qc": ["doublet_scores"]}, unit_rule="cluster_identity == 1")
    _, _, merge = commit_annotation_decisions(adata, [_decision(), technical], bundle, run_root=tmp_path, output_dir=tmp_path / "t")
    cells = adata.obs.cluster_identity.eq("1")
    assert set(adata.obs.loc[cells, "annotation_status"]) == {"excluded_technical"}
    assert set(adata.obs.loc[cells, "broad_annotation_status"]) == set(adata.obs.loc[cells, "fine_annotation_status"]) == {"not_applicable"}
    assert set(adata.obs.loc[cells, "annotation_technical_kind"]) == {"doublet"}
    review = adata.uns["annotation_review"]
    assert (review["technical_excluded"], review["reviewed_unresolved"], review["accepted"]) == (3, 0, 3)
    assert merge.action.tolist() == ["commit", "technical_exclusion"]
    assert merge.unit_rule.tolist()[1] == "cluster_identity == 1"
    assert review_coverage(adata)["technical_excluded"] == 3
    unbacked = validate_annotation_decision(_decision(cluster_id="1", disposition="technical_exclusion", technical_kind="ambient",
                                                      evidence_refs={"qc": ["manual:looks_ambient"]}), bundle)
    assert not unbacked.accepted
    with pytest.raises(ValueError, match="technical_kind"):
        _decision(disposition="technical_exclusion")


def test_fine_unresolved_notes_and_review_rows_are_recorded(tmp_path):
    adata = _adata(tmp_path)
    bundle = _bundle(adata, min_cells_per_sample=1)
    genus = _decision(final_label="T_Unresolved", identity_scope="fine_unresolved", confidence="high",
                      notes=["confidence_distribution high:2;low:1"])
    pending = _decision(cluster_id="1", evidence_refs={"deg": ["parent_markers::7"]}, review_required=True,
                        parent_unit="parent_cluster:7", unit_rule="split from 7")
    _, validation, merge = commit_annotation_decisions(adata, [genus, pending], bundle, run_root=tmp_path,
                                                       output_dir=tmp_path / "t", min_cells_high_confidence=1,
                                                       min_samples_high_confidence=1)
    assert validation.forced_review.tolist() == [False, False]
    first = adata.obs.cluster_identity.eq("0")
    assert set(adata.obs.loc[first, "broad_annotation_status"]) == {"accepted"}
    assert set(adata.obs.loc[first, "fine_annotation_status"]) == {"unresolved"}
    assert merge.action.tolist() == ["commit", "review_required"]
    assert merge.parent_unit.tolist() == [None, "parent_cluster:7"]
    assert "malignancy_status" not in adata.obs


def test_cell_layer_only_lowers_confidence(tmp_path):
    adata = _adata(tmp_path)
    bundle = _bundle(adata)
    layer = pd.DataFrame({"obs_name": ["c0", "c1", "c2", "c3"], "cell_confidence": ["low", "high", None, "low"],
                          "cell_state": ["Cycling", "", None, "Cycling"], "identity_source": ["knn_vote", None, None, "direct"]})
    pending = _decision(cluster_id="1", evidence_refs={"deg": ["parent_markers::7"]}, review_required=True)
    commit_annotation_decisions(adata, [_decision(), pending], bundle, run_root=tmp_path, output_dir=tmp_path / "t",
                                cell_layer=layer)
    assert adata.obs.annotation_confidence.tolist()[:4] == ["low", "medium", "medium", "unassigned"]
    assert adata.obs.cell_state.tolist()[:4] == ["Cycling", "", "", ""]
    assert adata.obs.annotation_identity_source.tolist()[:4] == ["knn_vote", "unit", "unit", ""]
    with pytest.raises(ValueError, match="absent from the input"):
        commit_annotation_decisions(adata, [_decision()], bundle, run_root=tmp_path, output_dir=tmp_path / "t",
                                    cell_layer=pd.DataFrame({"obs_name": ["zz"], "cell_confidence": ["low"]}))


def _run(tmp_path, monkeypatch, adata):
    run = tmp_path / "runs" / "test"
    (run / "config").mkdir(parents=True)
    (run / "artifacts").mkdir()
    config = run / "config" / "run.yaml"
    config.write_text(yaml.safe_dump({"run": {"run_id": "run1"}, "annotation": {"decision": {
        "sample_key": "sample_key", "min_cells_per_sample": 1, "min_cells_high_confidence": 3,
        "min_samples_high_confidence": 2}}}))
    adata.write_h5ad(run / "artifacts" / "adata_annotation_evidence.h5ad")
    monkeypatch.chdir(tmp_path)

    def call(module, *args):
        monkeypatch.setattr(sys, "argv", ["scoop", "--config", str(config), *args])
        return module.main()
    return run, call


def test_review_and_commit_bind_bundle_and_report_legacy(tmp_path, monkeypatch):
    run, call = _run(tmp_path, monkeypatch, _adata(tmp_path))
    out = run / "annotation" / "decisions"
    call(annotation_review)
    result = call(annotation_review)
    assert len(list((out / "history" / "annotation_evidence_bundle").glob("*.json"))) == 1
    assert result["coverage"][1]["processed"] == 1
    template = json.loads((out / "annotation_decision_template.json").read_text())["decisions"]
    assert [t["evidence_refs"]["deg"] for t in template] == [["cluster_markers"], []]
    assert template[0]["positive_markers"] == ["CD3E", "CD2", "LYZ"]
    decision = {**template[0], **_decision(confidence="high").to_dict(), "bundle_sha256": template[0]["bundle_sha256"]}
    decisions = tmp_path / "decisions.json"
    decisions.write_text(json.dumps([decision]))
    receipt = call(annotation_commit, "--decisions", str(decisions), "--evidence-bundle", str(out / "annotation_evidence_bundle.json"))
    assert {c["name"]: c["status"] for c in receipt["validation"]["checks"]}["bundle_partition_binding"] == "passed"
    committed = read_h5ad(run / "artifacts" / "adata_annotation_committed.h5ad")
    assert committed.obs.annotation_status.tolist()[:3] == ["accepted"] * 3

    legacy = json.loads((out / "annotation_evidence_bundle.json").read_text())
    for key in ("partition_sha256", "deg_tables", "composition_keys"):
        legacy.pop(key)
    (tmp_path / "legacy.json").write_text(json.dumps(legacy))
    stale = tmp_path / "stale.json"
    stale.write_text(json.dumps([decision]))
    with pytest.raises(ValueError, match="different evidence bundle"):
        call(annotation_commit, "--decisions", str(stale), "--evidence-bundle", str(tmp_path / "legacy.json"))
    decisions.write_text(json.dumps([{**decision, "bundle_sha256": None}]))
    receipt = call(annotation_commit, "--decisions", str(decisions), "--evidence-bundle", str(tmp_path / "legacy.json"))
    assert receipt["validation"]["status"] == "limited"
    assert {c["name"]: c["status"] for c in receipt["validation"]["checks"]}["bundle_partition_binding"] == "limited"
    assert read_h5ad(run / "artifacts" / "adata_annotation_committed.h5ad").obs.annotation_status.tolist()[:3] == ["accepted"] * 3


def test_unassessed_cnv_is_not_negative_and_cannot_be_a_cnv_basis(tmp_path):
    adata = _adata(tmp_path)
    adata.obs["fastcnv_tumor_evidence"] = ["cnv_unassessed"] * 3 + ["cnv_tumor_like", "cnv_normal_like", "cnv_unassessed"]
    adata.obs["fastcnv_threshold_source"] = "unassessed_insufficient_reference"
    bundle = _bundle(adata)
    unassessed, mixed = bundle.cluster_by_id()["0"].cnv_summary, bundle.cluster_by_id()["1"].cnv_summary
    assert (unassessed["has_tumor_evidence"], unassessed["fastcnv_assessed_cells"]) == (None, 0)
    assert (mixed["has_tumor_evidence"], mixed["fastcnv_unassessed_cells"]) == (True, 1)
    assert "fastcnv_threshold_source" in bundle.cluster_by_id()["0"].cnv_refs
    result = validate_annotation_decision(_decision(malignancy_status="malignant", malignancy_basis=["lineage_markers", "cnv"],
                                                    evidence_refs={"deg": ["cluster_markers"], "cnv": ["cnv_summary"]}), bundle)
    assert not result.accepted and any("not assessed" in r for r in result.reasons)


def test_standard_evidence_then_review_prefills_unit_deg(tmp_path, monkeypatch):
    from scoop.commands import annotation_evidence

    matrix = np.log1p(np.array([[6, 5, 4, 0, 0, 0]] * 4 + [[0, 0, 0, 4, 5, 6]] * 4, dtype=float))
    adata = AnnData(matrix, obs=pd.DataFrame({"cluster_identity": pd.Categorical(["0"] * 4 + ["1"] * 4),
                                              "sample_key": ["s1", "s2"] * 4}, index=[f"c{i}" for i in range(8)]),
                    var=pd.DataFrame(index=[f"g{i}" for i in range(6)]))
    adata.layers["log1p_norm"] = matrix.copy()
    run, call = _run(tmp_path, monkeypatch, adata)
    adata.write_h5ad(run / "artifacts" / "adata_programs.h5ad")
    call(annotation_evidence)
    call(annotation_review)
    out = run / "annotation" / "decisions"
    bundle = json.loads((out / "annotation_evidence_bundle.json").read_text())
    assert bundle["deg_tables"]["cluster_markers"]["label_key"] == "cluster_identity"
    template = json.loads((out / "annotation_decision_template.json").read_text())["decisions"]
    assert [t["evidence_refs"]["deg"] for t in template] == [["cluster_markers"], ["cluster_markers"]]
    # Markers come from the cited table; the three up genes per group tie on score.
    assert [set(t["positive_markers"][:3]) for t in template] == [{"g0", "g1", "g2"}, {"g3", "g4", "g5"}]


def test_commit_finds_the_cited_bundle_by_content(tmp_path, monkeypatch):
    run, call = _run(tmp_path, monkeypatch, _adata(tmp_path))
    out = run / "annotation" / "decisions"
    call(annotation_review)
    template = json.loads((out / "annotation_decision_template.json").read_text())["decisions"]
    decision = {**template[0], **_decision().to_dict(), "bundle_sha256": template[0]["bundle_sha256"]}
    decisions = tmp_path / "project" / "decisions.json"  # kept outside the review folder, as projects do
    decisions.parent.mkdir()
    decisions.write_text(json.dumps([decision]))
    evidence = read_h5ad(run / "artifacts" / "adata_annotation_evidence.h5ad")
    (tmp_path / "late.json").write_text("{}")
    register_file(evidence, key="late_report", path=tmp_path / "late.json", category="reports", schema="note.v1")
    evidence.write_h5ad(run / "artifacts" / "adata_annotation_evidence.h5ad")
    call(annotation_review)  # a changed bundle; the first one is archived byte-for-byte
    receipt = call(annotation_commit, "--decisions", str(decisions))
    checks = {c["name"]: c["status"] for c in receipt["validation"]["checks"]}
    assert checks["bundle_partition_binding"] == "passed"
    assert any("/history/annotation_evidence_bundle/" in path for path in receipt["inputs"])

    grown = read_h5ad(run / "artifacts" / "adata_annotation_evidence.h5ad")
    grown.obs["cluster_identity"] = ["0"] * 5 + ["1"]
    grown.write_h5ad(run / "artifacts" / "adata_annotation_evidence.h5ad")
    with pytest.raises(ValueError, match="different partition"):
        call(annotation_commit, "--decisions", str(decisions))
    with pytest.raises(FileNotFoundError):
        call(annotation_commit, "--decisions", str(decisions), "--evidence-bundle", str(tmp_path / "typo.json"))
    for path in out.rglob("annotation_evidence_bundle*.json"):
        path.unlink()
    for path in (out / "history" / "annotation_evidence_bundle").glob("*.json"):
        path.unlink()
    with pytest.raises(ValueError, match="different evidence bundle"):
        call(annotation_commit, "--decisions", str(decisions))
    decisions.write_text(json.dumps([{**decision, "bundle_sha256": None}]))
    receipt = call(annotation_commit, "--decisions", str(decisions))
    assert {c["name"]: c["status"] for c in receipt["validation"]["checks"]}["bundle_partition_binding"] == "not_assessed"
    assert receipt["validation"]["status"] == "limited"


def test_group_citations_bind_to_a_declared_partition(tmp_path):
    adata = _adata(tmp_path)
    pd.DataFrame({"group": ["0", "1"], "names": ["CD3E", "LYZ"], "scores": [2.0, 1.0]}).to_csv(
        tmp_path / "unkeyed.tsv", sep="\t", index=False)
    register_file(adata, key="unkeyed_markers", path=tmp_path / "unkeyed.tsv", schema="cluster_markers.v1")
    bundle = _bundle(adata)
    sibling = validate_annotation_decision(_decision(cluster_id="1", evidence_refs={"deg": ["cluster_markers::0"]}), bundle)
    assert not sibling.accepted  # a sibling unit of the same partition is other cells
    assert validate_annotation_decision(_decision(evidence_refs={"deg": ["cluster_markers::0"]}), bundle).accepted
    unkeyed = validate_annotation_decision(_decision(evidence_refs={"deg": ["unkeyed_markers::0"]}), bundle)
    assert not unkeyed.accepted and any("declares no label_key" in r for r in unkeyed.reasons)


def test_tsv_groups_stay_labels(tmp_path):
    adata = AnnData(np.ones((2, 1)), obs=pd.DataFrame({"cluster_identity": ["01", "02"]}, index=["a", "b"]))
    pd.DataFrame({"group": ["01", "02"], "names": ["G", "H"], "scores": ["3", "1"]}).to_csv(
        tmp_path / "deg.tsv", sep="\t", index=False)
    register_file(adata, key="cluster_markers", path=tmp_path / "deg.tsv", schema="cluster_markers.v1",
                  extra={"label_key": "cluster_identity"})
    adata.uns["top_markers_01"] = ["FROM_ANOTHER_PARTITION"]
    clusters = _bundle(adata).cluster_by_id()
    assert clusters["01"].deg_refs == ("cluster_markers",) and clusters["01"].top_markers == ("G",)


def test_malignancy_gate_closes_status_bypasses(tmp_path):
    bundle = _bundle(_adata(tmp_path), min_cells_per_sample=1)
    malignant = dict(parent_label="Malignant", canonical_label="GBM_Malignant", final_label="GBM_Malignant_AC")
    assert not validate_annotation_decision(_decision(**malignant, malignancy_status="not_applicable"), bundle).accepted
    high = validate_annotation_decision(_decision(malignancy_status="ambiguous", confidence="high"), bundle,
                                        min_cells_high_confidence=1)
    assert not high.accepted and any("ambiguous" in r for r in high.reasons)
    assert validate_annotation_decision(_decision(malignancy_status="ambiguous"), bundle).accepted


def test_recommit_clears_claims_of_a_superseded_decision(tmp_path):
    adata = _adata(tmp_path)
    bundle = _bundle(adata)
    malignant = _decision(final_label="GBM_Malignant", malignancy_status="malignant", malignancy_basis=["lineage_markers"])
    commit_annotation_decisions(adata, [malignant], bundle, run_root=tmp_path, output_dir=tmp_path / "t")
    assert set(adata.obs.malignancy_status[:3]) == {"malignant"}
    commit_annotation_decisions(adata, [_decision(review_required=True)], bundle, run_root=tmp_path, output_dir=tmp_path / "t")
    assert set(adata.obs.annotation_status[:3]) == {"review_required"}
    assert set(adata.obs.malignancy_status[:3]) == {"not_assessed"} and set(adata.obs.malignancy_basis[:3]) == {""}


def test_technical_exclusion_needs_artefact_evidence_not_sample_support(tmp_path):
    adata = _adata(tmp_path)
    (tmp_path / "template.tsv").write_text("x\n")
    register_file(adata, key="annotation_evidence_template", path=tmp_path / "template.tsv", schema="annotation_evidence.v1")
    bundle = _bundle(adata, min_cells_per_sample=1)
    fields = dict(cluster_id="1", disposition="technical_exclusion", technical_kind="doublet", confidence="high")
    blank = validate_annotation_decision(_decision(**fields, evidence_refs={"qc": ["annotation_evidence_template"]}), bundle,
                                         min_cells_high_confidence=1, min_samples_high_confidence=5)
    assert not blank.accepted
    scored = validate_annotation_decision(_decision(**fields, evidence_refs={"qc": ["doublet_scores"]}), bundle,
                                          min_cells_high_confidence=1, min_samples_high_confidence=5)
    assert scored.accepted, scored.reasons


def test_decisions_from_tables_parse_nulls_and_joined_lists():
    decision = _decision(identity_scope=None, disposition=None, malignancy_status="malignant",
                         malignancy_basis="lineage_markers;cnv", evidence_refs={"deg": "cluster_markers"})
    assert (decision.identity_scope, decision.disposition) == ("broad", "identity")
    assert decision.malignancy_basis == ("lineage_markers", "cnv")
    assert decision.evidence_refs["deg"] == ("cluster_markers",)


def test_marker_deg_refuses_duplicate_gene_names():
    from scsp_agent_sop.annotation import run_markers_scanpy

    adata = AnnData(np.ones((4, 2)), obs=pd.DataFrame({"g": ["a", "a", "b", "b"]}, index=list("wxyz")),
                    var=pd.DataFrame(index=["SAME", "SAME"]))
    with pytest.raises(ValueError, match="unique var_names"):
        run_markers_scanpy(adata, groupby="g", layer="X")
