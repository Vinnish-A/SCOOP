from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any


@dataclass(frozen=True)
class EngineSpec:
    engine_id: str
    task_type: str
    version: str
    input_schema: str
    output_schema: str
    consumes: tuple[str, ...]
    produces: tuple[str, ...]
    default_cli: tuple[str, ...]
    writes_h5ad_fields: tuple[str, ...]
    writes_external_artifacts: tuple[str, ...]
    quality_gates: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "EngineSpec":
        tuple_fields = {
            "consumes",
            "produces",
            "default_cli",
            "writes_h5ad_fields",
            "writes_external_artifacts",
            "quality_gates",
        }
        normalized = dict(data)
        for field in tuple_fields:
            normalized[field] = tuple(normalized.get(field, ()))
        return cls(**normalized)
