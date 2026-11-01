from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping

from .decision_schema import ALLOWED_CONFIDENCE, AnnotationDecision
from .evidence_bundle import AnnotationEvidenceBundle, ClusterEvidence


TUMOR_TERMS = ("malignant", "tumor", "tumour", "cancer", "neoplastic")
# Defaults for the tumor_gate block of configs/annotation_decision_schema.yaml.
# Lineage markers decide malignancy; CNV is supporting evidence, never required.
TUMOR_GATE = {
    "enabled": True,
    "tumor_terms": TUMOR_TERMS,
    "malignancy_status": ("malignant", "probable_malignant", "normal", "ambiguous", "not_applicable"),
    "malignant_status": ("malignant", "probable_malignant"),
    "malignancy_basis": ("lineage_markers", "tumour_programmes", "cnv", "snv", "atac", "reference_mapping", "histology"),
    "required_basis": ("lineage_markers",),
    "cnv_schema_prefixes": ("cnv", "fastcnvpy"),
}
DECLARED_EXTERNAL_PREFIXES = ("external:", "manual:")
# Registered schemas that can evidence a technical exclusion besides cluster DEG.
TECHNICAL_EVIDENCE_PREFIXES = ("qc", "scrublet", "doublet", "ambient", "core_quality")
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
    tumor_gate: Mapping[str, Any] | None = None,
) -> ValidationResult:
    del skills
    gate = {**TUMOR_GATE, **dict(tumor_gate or {})}
    technical = decision.disposition == "technical_exclusion"
    reasons: list[str] = []
    forced_review = False
    if decision.identity_scope == "unresolved":
        reasons.append("identity explicitly unresolved")
        forced_review = True
    if decision.cluster_key != evidence_bundle.cluster_key or decision.run_id != evidence_bundle.run_id:
        return ValidationResult(decision.cluster_id, False, False, True, decision.confidence,
                                ("decision run/cluster axis differs from evidence bundle",), True)
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
    valid_deg, uncovered_deg = _deg_citations(decision, evidence_bundle)
    if uncovered_deg:
        reasons.append(f"cited DEG does not test this unit (no such label_key/group): {', '.join(uncovered_deg)}")
        forced_review = True
    label_text = " ".join(
        str(value or "") for value in (decision.parent_label, decision.canonical_label, decision.final_label)
    ).lower()
    if technical:
        # An exclusion makes no identity claim, but it must cite evidence of the
        # artefact: unit DEG or a registered QC/doublet table, not any file.
        if not valid_deg and not _cites_registered_schema(decision, evidence_bundle, TECHNICAL_EVIDENCE_PREFIXES):
            reasons.append("technical exclusion requires unit DEG or a cited QC/doublet artefact registered with such a schema")
            forced_review = True
    else:
        if not valid_deg:
            reasons.append("cluster differential-gene evidence is not cited; marker scores or transferred labels alone cannot qualify identity")
            forced_review = True
        state_text = " ".join(str(value or "") for value in (decision.cell_state, decision.functional_modifier)).lower()
        if any(term in label_text for term in PROGRAM_ONLY_TERMS) and not any(term in state_text for term in PROGRAM_ONLY_TERMS):
            reasons.append("program-dominated label should be represented as cell_state or functional_modifier")
            forced_review = True
    if gate["enabled"]:
        gate_reasons = _malignancy_gate(decision, evidence_bundle, cluster, gate, label_text, bool(valid_deg))
        reasons.extend(gate_reasons)
        forced_review = forced_review or bool(gate_reasons)
    if decision.conflicts and decision.confidence == "high":
        reasons.append("marker conflicts prevent high confidence")
        forced_review = True
    if decision.confidence == "high":
        if cluster.n_cells < min_cells_high_confidence:
            reasons.append("cluster has too few cells for high confidence")
            forced_review = True
        # Sample support is an identity criterion; a single-library artefact is still an artefact.
        if not technical and cluster.n_samples is None:
            reasons.append("sample support is unknown (sample key absent); high confidence requires it")
            forced_review = True
        elif not technical and cluster.n_samples < min_samples_high_confidence:
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


def _deg_citations(decision: AnnotationDecision, bundle: AnnotationEvidenceBundle) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Resolve DEG citations to (table, group) pairs the cited table really tests.

    "ref" means group == cluster_id in a table whose label_key is the bundle's
    cluster_key; "ref::group" names a group of another partition (or this
    unit itself). A table that declares no label_key binds nothing. Declared
    external:/manual: refs are the analyst's own, unbound DEG. A "markers"
    citation counts only when it resolves; an unresolved "deg" one is flagged.
    """
    valid, uncovered = [], []
    for category in ("deg", "markers"):
        for ref in map(str, decision.evidence_refs.get(category, ())):
            if ref.startswith(DECLARED_EXTERNAL_PREFIXES):
                if category == "deg":
                    valid.append(ref)
                continue
            name, qualified, group = ref.partition("::")
            table = bundle.deg_tables.get(name)
            if table is None and category == "markers":
                continue
            if table is not None and table["label_key"] is None:
                # Without a declared partition the group cannot be bound to anything.
                uncovered.append(f"{ref} (table declares no label_key; re-register it with the partition it tests)")
                continue
            own_partition = table is not None and table["label_key"] == bundle.cluster_key
            if qualified:
                # In the unit's own partition "::group" may name only the unit
                # itself: a sibling group is a different set of cells.
                covered = (table is not None and group in table["groups"]
                           and (not own_partition or group == str(decision.cluster_id)))
            else:
                covered = own_partition and str(decision.cluster_id) in table["groups"]
            (valid if covered else uncovered).append(ref)
    return tuple(valid), tuple(sorted(set(uncovered)))


def _registered_refs(decision: AnnotationDecision, bundle: AnnotationEvidenceBundle) -> tuple[str, ...]:
    registered = {str(key) for category in ("tables", "artifacts", "reports", "figures")
                  for key in bundle.file_registry.get(category, {})}
    return tuple(ref for refs in decision.evidence_refs.values() for ref in map(str, refs)
                 if ref.split("::", 1)[0] in registered)


def _malignancy_gate(decision: AnnotationDecision, bundle: AnnotationEvidenceBundle, cluster: ClusterEvidence,
                     gate: Mapping[str, Any], label_text: str, has_valid_deg: bool) -> list[str]:
    """Malignancy is declared, not read from label text, and rests on lineage markers."""
    reasons = []
    status, basis = decision.malignancy_status, set(decision.malignancy_basis)
    tumor_label = any(term in label_text for term in gate["tumor_terms"])
    if status is not None and status not in gate["malignancy_status"]:
        reasons.append(f"unknown malignancy_status {status!r}")
    unknown_basis = sorted(basis - set(gate["malignancy_basis"]))
    if unknown_basis:
        reasons.append(f"unknown malignancy_basis: {', '.join(unknown_basis)}")
    if status is None and tumor_label:
        reasons.append("tumour term in label requires a declared malignancy_status")
    if status in {"normal", "not_applicable"} and tumor_label:
        reasons.append(f"malignancy_status {status} contradicts the tumour term in the label")
    if status == "ambiguous" and decision.confidence == "high":
        reasons.append("ambiguous malignancy cannot be high confidence")
    if status in gate["malignant_status"]:
        missing = sorted(set(gate["required_basis"]) - basis)
        if missing or not has_valid_deg:
            reasons.append(f"malignancy_status {status} requires basis {', '.join(gate['required_basis'])} "
                           "and a valid DEG citation; CNV is supporting only")
    if "cnv" in basis and not _cites_registered_schema(decision, bundle, tuple(gate["cnv_schema_prefixes"])):
        reasons.append("malignancy_basis cnv requires a cited artefact registered with a CNV schema")
    if "cnv" in basis and cluster.cnv_summary.get("fastcnv_assessed_cells") == 0:
        reasons.append("malignancy_basis cnv cites CNV that was not assessed for this unit (cnv_unassessed)")
    return reasons


def _cites_registered_schema(decision: AnnotationDecision, bundle: AnnotationEvidenceBundle, prefixes: tuple[str, ...]) -> bool:
    records = {str(key): record for category in bundle.file_registry.values() if isinstance(category, Mapping)
               for key, record in category.items() if isinstance(record, Mapping)}
    return any(str(records.get(ref, {}).get("schema", "")).startswith(prefixes) for ref in _registered_refs(decision, bundle))


def _unknown_evidence_refs(decision: AnnotationDecision, bundle: AnnotationEvidenceBundle, cluster: ClusterEvidence) -> tuple[str, ...]:
    known = set(cluster.marker_refs) | set(cluster.nmf_refs) | set(cluster.reference_refs) | set(cluster.cnv_refs)
    for category in ("tables", "artifacts", "reports", "figures"):
        known.update(map(str, bundle.file_registry.get(category, {}).keys()))
    unknown = []
    for refs in decision.evidence_refs.values():
        for ref in refs:
            if str(ref).split("::", 1)[0] not in known and not str(ref).startswith(("external:", "manual:", "skill:")):
                unknown.append(str(ref))
    return tuple(sorted(set(unknown)))
