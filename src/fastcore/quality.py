from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np


@dataclass(frozen=True)
class CoreQualitySummary:
    accepted: bool
    metrics: dict[str, float] = field(default_factory=dict)
    reasons: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def pc_subspace_cosine(x: np.ndarray, y: np.ndarray, n_pcs: int | None = None) -> float:
    k = min(x.shape[1], y.shape[1], n_pcs or min(x.shape[1], y.shape[1]))
    if k <= 0:
        return 0.0
    qx, _ = np.linalg.qr(np.asarray(x[:, :k], dtype=float))
    qy, _ = np.linalg.qr(np.asarray(y[:, :k], dtype=float))
    singular = np.linalg.svd(qx.T @ qy, compute_uv=False)
    return float(np.mean(np.clip(singular, 0.0, 1.0)))


def graph_density(graph) -> float:
    n = int(graph.shape[0])
    if n <= 1:
        return 0.0
    nnz = int(graph.nnz) if hasattr(graph, "nnz") else int(np.count_nonzero(graph))
    return float(nnz / (n * (n - 1)))
