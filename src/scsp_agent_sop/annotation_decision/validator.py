from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping

from .decision_schema import ALLOWED_CONFIDENCE, AnnotationDecision
from .evidence_bundle import AnnotationEvidenceBundle, ClusterEvidence


TUMOR_TERMS = ("malignant", "tumor", "tumour", "cancer", "neoplastic")
PROGRAM_ONLY_TERMS = ("ribosomal", "mitochondrial", "stress", "cell-cycle", "cell_cycle", "cycling")


@dataclass(frozen=True)
class ValidationResult:
    cluster_id: str
    valid: bool
    accepted: bool
    review_required: bool
    confidence: str | None
    reasons: tuple[str, ...]
    forced_review: bool = False

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def validate_annotation_decision(
    decision: AnnotationDecision,
    evidence_bundle: AnnotationEvidenceBundle,
    skills: Mapping[str, Any] | None = None,
    *,
    min_cells_high_confidence: int = 20,
    min_samples_high_confidence: int = 2,
) -> ValidationResult:
    del skills
    reasons: list[str] = []
    forced_review = False
    cluster = evidence_bundle.cluster_by_id().get(str(decision.cluster_id))
    if cluster is None:
        return ValidationResult(decision.cluster_id, False, False, True, decision.confidence, ("cluster not found in evidence bundle",), True)
    if decision.confidence not in ALLOWED_CONFIDENCE:
        return ValidationResult(decision.cluster_id, False, False, True, decision.confidence, ("invalid confidence",), True)
    if not decision.final_label.strip():
        return ValidationResult(decision.cluster_id, False, False, True, decision.confidence, ("empty final_label",), True)
    unknown_refs = _unknown_evidence_refs(decision, evidence_bundle, cluster)
    if unknown_refs:
        reasons.append(f"unknown evidence refs: {', '.join(unknown_refs)}")
        forced_review = True
    label_text = " ".join(
        str(value or "") for value in (decision.parent_label, decision.canonical_label, decision.final_label)
    ).lower()
    state_text = " ".join(str(value or "") for value in (decision.cell_state, decision.functional_modifier)).lower()
    if any(term in label_text for term in PROGRAM_ONLY_TERMS) and not any(term in state_text for term in PROGRAM_ONLY_TERMS):
        reasons.append("program-dominated label should be represented as cell_state or functional_modifier")
        forced_review = True
    if evidence_bundle.is_tumor and _uses_tumor_label(decision):
        has_cnv = _decision_cites_cnv(decision) or bool(cluster.cnv_refs and cluster.cnv_summary.get("has_tumor_evidence"))
        if not has_cnv:
            reasons.append("tumor/malignant label requires CNV evidence or review")
            forced_review = True
    if decision.conflicts and decision.confidence == "high":
        reasons.append("marker conflicts prevent high confidence")
        forced_review = True
    if decision.confidence == "high":
        if cluster.n_cells < min_cells_high_confidence:
            reasons.append("cluster has too few cells for high confidence")
            forced_review = True
        if cluster.n_samples is not None and cluster.n_samples < min_samples_high_confidence:
            reasons.append("cluster has too few samples for high confidence")
            forced_review = True
    review_required = bool(decision.review_required or forced_review)
    accepted = not forced_review and not review_required
    if not reasons and review_required:
        reasons.append("decision requested review")
    if not reasons:
        reasons.append("validation passed")
    return ValidationResult(
        cluster_id=decision.cluster_id,
        valid=True,
        accepted=accepted,
        review_required=review_required,
        confidence=decision.confidence,
        reasons=tuple(reasons),
        forced_review=forced_review,
    )


def _uses_tumor_label(decision: AnnotationDecision) -> bool:
    text = " ".join(str(value or "") for value in (decision.parent_label, decision.canonical_label, decision.final_label)).lower()
    return any(term in text for term in TUMOR_TERMS)


def _decision_cites_cnv(decision: AnnotationDecision) -> bool:
    for key, refs in decision.evidence_refs.items():
        if "cnv" in key.lower() and refs:
            return True
        if any("cnv" in str(ref).lower() or "fastcnv" in str(ref).lower() for ref in refs):
            return True
    return False


def _unknown_evidence_refs(decision: AnnotationDecision, bundle: AnnotationEvidenceBundle, cluster: ClusterEvidence) -> tuple[str, ...]:
    known = set(cluster.marker_refs) | set(cluster.nmf_refs) | set(cluster.reference_refs) | set(cluster.cnv_refs)
    for category in ("tables", "artifacts", "reports", "figures"):
        known.update(map(str, bundle.file_registry.get(category, {}).keys()))
    unknown = []
    for refs in decision.evidence_refs.values():
        for ref in refs:
            if ref not in known and not str(ref).startswith(("external:", "manual:", "skill:")):
                unknown.append(str(ref))
    return tuple(sorted(set(unknown)))
