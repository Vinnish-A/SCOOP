from __future__ import annotations

from .backend_plan import FastCorePlan, plan_fastcore_backend
from .runtime import FastCoreCapabilities, detect_capabilities

__all__ = [
    "FastCoreCapabilities",
    "FastCorePlan",
    "detect_capabilities",
    "plan_fastcore_backend",
]
