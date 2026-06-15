from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any
import json


@dataclass(frozen=True)
class ArtifactBundle:
    schema_version: str
    engine_id: str
    task_type: str
    run_id: str
    status: str
    inputs: dict[str, Any]
    outputs: dict[str, Any]
    quality: dict[str, Any]
    timings: dict[str, float]
    registry_patch: dict[str, Any]
    decision_log_patch: dict[str, Any]
    review_required: bool = False
    review_reason: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ArtifactBundle":
        normalized = dict(data)
        normalized["timings"] = {str(k): float(v) for k, v in normalized.get("timings", {}).items()}
        normalized["review_required"] = bool(normalized.get("review_required", False))
        return cls(**normalized)

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), ensure_ascii=False, sort_keys=True)

    @classmethod
    def from_json(cls, text: str) -> "ArtifactBundle":
        return cls.from_dict(json.loads(text))

    def write_json(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return path

    @classmethod
    def read_json(cls, path: str | Path) -> "ArtifactBundle":
        return cls.from_json(Path(path).read_text(encoding="utf-8"))
