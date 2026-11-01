from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping
import json


ALLOWED_CONFIDENCE = ("high", "medium", "low")
ALLOWED_IDENTITY_SCOPE = ("broad", "fine", "unresolved", "fine_unresolved")
ALLOWED_DISPOSITION = ("identity", "technical_exclusion")
ALLOWED_TECHNICAL_KIND = ("doublet", "multiplet", "ambient", "low_quality", "artefact")


@dataclass(frozen=True)
class AnnotationDecision:
    schema_version: str
    run_id: str
    cluster_id: str
    cluster_key: str
    parent_label: str
    canonical_label: str
    cell_state: str | None
    functional_modifier: str | None
    final_label: str
    confidence: str
    evidence_refs: dict[str, tuple[str, ...]]
    positive_markers: tuple[str, ...]
    negative_markers_absent: tuple[str, ...]
    conflicts: tuple[str, ...]
    review_required: bool
    reason: str
    identity_scope: str = "broad"
    # fine_unresolved: broad lineage accepted, fine identity reviewed but unresolved.
    disposition: str = "identity"
    technical_kind: str | None = None
    # Never gates; conflicts stays reserved for genuine marker conflicts.
    notes: tuple[str, ...] = ()
    # Declared, not inferred from label text; vocabularies live in tumor_gate.
    malignancy_status: str | None = None
    malignancy_basis: tuple[str, ...] = ()
    parent_unit: str | None = None
    unit_rule: str | None = None
    bundle_sha256: str | None = None

    def __post_init__(self) -> None:
        if self.confidence not in ALLOWED_CONFIDENCE:
            raise ValueError(f"confidence must be one of {ALLOWED_CONFIDENCE}: {self.confidence!r}")
        if not str(self.final_label).strip():
            raise ValueError("final_label must not be empty")
        if self.identity_scope not in ALLOWED_IDENTITY_SCOPE:
            raise ValueError(f"identity_scope must be one of {ALLOWED_IDENTITY_SCOPE}")
        if self.disposition not in ALLOWED_DISPOSITION:
            raise ValueError(f"disposition must be one of {ALLOWED_DISPOSITION}")
        if (self.disposition == "technical_exclusion") != (self.technical_kind is not None):
            raise ValueError("technical_kind is required for, and only for, disposition technical_exclusion")
        if self.technical_kind is not None and self.technical_kind not in ALLOWED_TECHNICAL_KIND:
            raise ValueError(f"technical_kind must be one of {ALLOWED_TECHNICAL_KIND}")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "AnnotationDecision":
        normalized = dict(data)
        normalized["evidence_refs"] = {
            str(key): tuple(map(str, values or ())) for key, values in dict(normalized.get("evidence_refs", {})).items()
        }
        for field in ("positive_markers", "negative_markers_absent", "conflicts", "notes", "malignancy_basis"):
            normalized[field] = tuple(map(str, normalized.get(field) or ()))
        normalized["review_required"] = bool(normalized.get("review_required", False))
        normalized.setdefault("schema_version", "scoop.annotation_decision.v1")
        normalized.setdefault("cell_state", None)
        normalized.setdefault("functional_modifier", None)
        return cls(**normalized)

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), ensure_ascii=False, sort_keys=True)

    @classmethod
    def from_json(cls, text: str) -> "AnnotationDecision":
        return cls.from_dict(json.loads(text))


def decisions_from_json(text: str) -> list[AnnotationDecision]:
    payload = json.loads(text)
    if isinstance(payload, dict) and "decisions" in payload:
        payload = payload["decisions"]
    if isinstance(payload, dict):
        payload = [payload]
    if not isinstance(payload, list):
        raise ValueError("decision JSON must be a decision object, a list, or {'decisions': [...]}")
    return [AnnotationDecision.from_dict(item) for item in payload]
