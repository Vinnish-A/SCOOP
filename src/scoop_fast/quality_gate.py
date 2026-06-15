from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any


@dataclass(frozen=True)
class QualityGateResult:
    gate_id: str
    passed: bool
    severity: str = "info"
    reason: str = ""
    metrics: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["metrics"] = dict(self.metrics or {})
        return data
