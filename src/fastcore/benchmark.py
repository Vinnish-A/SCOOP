from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass(frozen=True)
class CoreBenchmarkRecord:
    backend: str
    n_obs: int
    n_vars: int
    wall_time_s: float
    peak_rss_mb: float | None = None
    peak_gpu_mb: float | None = None
    quality: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
