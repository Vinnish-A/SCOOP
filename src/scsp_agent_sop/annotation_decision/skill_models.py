from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class AnnotationSkill:
    skill_id: str
    version: str
    species: str
    tissue: str | None = None
    disease_context: str | None = None
    lineage_scope: str | None = None
    gene_symbol_namespace: str = "HGNC"
    sources: tuple[str, ...] = ()
    naming_policy: dict[str, Any] | None = None

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "AnnotationSkill":
        normalized = dict(data)
        normalized["sources"] = tuple(map(str, normalized.get("sources", ())))
        normalized["naming_policy"] = dict(normalized.get("naming_policy", {}))
        return cls(**normalized)
