from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping

from scsp_agent_sop.config import deep_get

from .runtime import FastCoreCapabilities, detect_capabilities


FASTCORE_BACKENDS = (
    "omicverse_rust_oom",
    "omicverse_gpu_rapids",
    "omicverse_cpu_gpu_mixed",
    "omicverse_cpu",
)


@dataclass(frozen=True)
class FastCorePlan:
    selected_backend: str
    fallback_required: bool
    fallback_backend: str
    reasons: list[str] = field(default_factory=list)
    capabilities: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _input_is_backed(input_path: str | Path | None, adata: Any | None) -> bool:
    if getattr(adata, "isbacked", False):
        return True
    if input_path is None:
        return False
    return str(input_path).endswith(".h5ad")


def plan_fastcore_backend(
    cfg: Mapping[str, Any],
    *,
    adata: Any | None = None,
    input_path: str | Path | None = None,
    capabilities: FastCoreCapabilities | None = None,
) -> FastCorePlan:
    """Choose one backend before execution; do not chain per-step fallbacks."""
    caps = capabilities or detect_capabilities()
    fallback_backend = str(deep_get(cfg, "core.fastcore.fallback_backend", deep_get(cfg, "core.fallback_engine", "scanpy_legacy")))
    allowed = tuple(deep_get(cfg, "core.fastcore.allowed_backends", FASTCORE_BACKENDS))
    policy = str(deep_get(cfg, "core.fastcore.backend_policy", "auto"))
    reasons = list(caps.reasons)

    if policy != "auto":
        if policy in allowed and _backend_capable(policy, caps):
            return FastCorePlan(policy, False, fallback_backend, reasons, caps.to_dict())
        reasons.append(f"requested_backend_unavailable:{policy}")
        return FastCorePlan(fallback_backend, True, fallback_backend, reasons, caps.to_dict())

    has_cpu_backend = caps.omicverse_available or getattr(caps, "vendored_omicverse_available", False)
    if not has_cpu_backend:
        return FastCorePlan(fallback_backend, True, fallback_backend, reasons, caps.to_dict())
    if not bool(deep_get(cfg, "core.fastcore.enable_omicverse_cpu_backend", True)):
        reasons.append("omicverse_cpu_backend_disabled")
        return FastCorePlan(fallback_backend, True, fallback_backend, reasons, caps.to_dict())

    n_obs = int(getattr(adata, "n_obs", 0) or 0)
    n_vars = int(getattr(adata, "n_vars", 0) or 0)
    nnz = int(getattr(getattr(adata, "X", None), "nnz", 0) or 0)
    large_cells = n_obs >= int(deep_get(cfg, "core.fastcore.auto.rust_oom_min_cells", 100000))
    large_nnz = nnz >= int(deep_get(cfg, "core.fastcore.auto.rust_oom_min_nnz", 300000000))
    gpu_cells = n_obs >= int(deep_get(cfg, "core.fastcore.auto.gpu_min_cells", 30000))
    backed = _input_is_backed(input_path, adata)

    if (
        "omicverse_rust_oom" in allowed
        and caps.rust_backend_available
        and (backed or large_cells or large_nnz)
        and bool(deep_get(cfg, "core.fastcore.auto.prefer_rust_when_backed", True))
    ):
        return FastCorePlan("omicverse_rust_oom", False, fallback_backend, reasons, caps.to_dict())
    if "omicverse_gpu_rapids" in allowed and caps.cuda_available and caps.rapids_available and gpu_cells:
        return FastCorePlan("omicverse_gpu_rapids", False, fallback_backend, reasons, caps.to_dict())
    if (
        "omicverse_cpu_gpu_mixed" in allowed
        and caps.cuda_available
        and caps.torch_available
        and bool(deep_get(cfg, "core.fastcore.auto.prefer_gpu_when_available", True))
    ):
        return FastCorePlan("omicverse_cpu_gpu_mixed", False, fallback_backend, reasons, caps.to_dict())
    if "omicverse_cpu" in allowed and has_cpu_backend:
        return FastCorePlan("omicverse_cpu", False, fallback_backend, reasons, caps.to_dict())
    reasons.append("no_allowed_omicverse_backend_available")
    return FastCorePlan(fallback_backend, True, fallback_backend, reasons, caps.to_dict())


def _backend_capable(backend: str, caps: FastCoreCapabilities) -> bool:
    if backend == "scanpy_legacy":
        return True
    if backend == "omicverse_cpu":
        return caps.omicverse_available or getattr(caps, "vendored_omicverse_available", False)
    if backend == "omicverse_cpu_gpu_mixed":
        return caps.omicverse_available and caps.torch_available and caps.cuda_available
    if backend == "omicverse_gpu_rapids":
        return caps.omicverse_available and caps.cuda_available and caps.rapids_available
    if backend == "omicverse_rust_oom":
        return caps.omicverse_available and caps.rust_backend_available
    return False
