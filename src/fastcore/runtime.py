from __future__ import annotations

from dataclasses import asdict, dataclass, field
from importlib import util
from typing import Any


def _module_available(name: str) -> bool:
    return util.find_spec(name) is not None


@dataclass(frozen=True)
class FastCoreCapabilities:
    omicverse_available: bool
    torch_available: bool
    cuda_available: bool
    rapids_available: bool
    anndataoom_available: bool
    rust_backend_available: bool
    selected_backend: str | None = None
    fallback_required: bool = False
    reasons: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def detect_capabilities() -> FastCoreCapabilities:
    """Detect optional FastCore runtime features without importing heavy stacks eagerly."""
    reasons: list[str] = []
    omicverse_available = _module_available("omicverse")
    torch_available = _module_available("torch")
    cuda_available = False
    if torch_available:
        try:
            import torch

            cuda_available = bool(torch.cuda.is_available())
        except Exception as exc:  # pragma: no cover - environment-dependent.
            reasons.append(f"torch_cuda_probe_failed:{type(exc).__name__}")
    rapids_available = any(_module_available(name) for name in ("rapids_singlecell", "cuml", "cugraph"))
    anndataoom_available = any(_module_available(name) for name in ("anndata_oom", "anndata_rs"))
    rust_backend_available = omicverse_available and anndataoom_available
    if not omicverse_available:
        reasons.append("omicverse_unavailable")
    if torch_available and not cuda_available:
        reasons.append("torch_without_cuda")
    if cuda_available and not rapids_available:
        reasons.append("rapids_unavailable")
    if not anndataoom_available:
        reasons.append("anndataoom_unavailable")
    return FastCoreCapabilities(
        omicverse_available=omicverse_available,
        torch_available=torch_available,
        cuda_available=cuda_available,
        rapids_available=rapids_available,
        anndataoom_available=anndataoom_available,
        rust_backend_available=rust_backend_available,
        reasons=reasons,
    )
