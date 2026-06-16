from __future__ import annotations

import numpy as np
from anndata import AnnData

from fastcore.backend_plan import FastCorePlan, plan_fastcore_backend
from fastcore.runtime import FastCoreCapabilities, detect_capabilities


def test_detect_capabilities_does_not_raise():
    caps = detect_capabilities()
    data = caps.to_dict()
    assert "omicverse_available" in data
    assert "reasons" in data


def test_planner_falls_back_when_vendored_backend_disabled():
    cfg = {"core": {"engine": "fastcore", "fastcore": {"fallback_backend": "scanpy_legacy", "enable_fastcore_cpu_backend": False}}}
    adata = AnnData(np.ones((4, 3)))
    caps = FastCoreCapabilities(
        omicverse_available=False,
        torch_available=False,
        cuda_available=False,
        anndataoom_available=False,
        rust_backend_available=False,
        vendored_omicverse_available=True,
        reasons=["omicverse_unavailable"],
    )
    plan = plan_fastcore_backend(cfg, adata=adata, capabilities=caps)
    assert plan == FastCorePlan(
        selected_backend="scanpy_legacy",
        fallback_required=True,
        fallback_backend="scanpy_legacy",
        reasons=["omicverse_unavailable", "fastcore_cpu_backend_disabled"],
        capabilities=caps.to_dict(),
    )


def test_planner_selects_vendored_cpu_backend_by_default():
    cfg = {"core": {"engine": "fastcore", "fastcore": {"fallback_backend": "scanpy_legacy"}}}
    adata = AnnData(np.ones((4, 3)))
    caps = FastCoreCapabilities(
        omicverse_available=True,
        torch_available=False,
        cuda_available=False,
        anndataoom_available=False,
        rust_backend_available=False,
        vendored_omicverse_available=True,
        reasons=[],
    )
    plan = plan_fastcore_backend(cfg, adata=adata, capabilities=caps)
    assert plan.selected_backend == "fastcore_cpu"
    assert plan.fallback_required is False


def test_planner_selects_cpu_when_enabled_and_available():
    cfg = {
        "core": {
            "engine": "fastcore",
            "fastcore": {
                "fallback_backend": "scanpy_legacy",
                "enable_fastcore_cpu_backend": True,
                "allowed_backends": ["fastcore_cpu"],
            },
        }
    }
    adata = AnnData(np.ones((4, 3)))
    caps = FastCoreCapabilities(
        omicverse_available=True,
        torch_available=False,
        cuda_available=False,
        anndataoom_available=False,
        rust_backend_available=False,
        vendored_omicverse_available=True,
        reasons=[],
    )
    plan = plan_fastcore_backend(cfg, adata=adata, capabilities=caps)
    assert plan.selected_backend == "fastcore_cpu"
    assert plan.fallback_required is False
