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


def test_planner_falls_back_when_omicverse_missing():
    cfg = {"core": {"engine": "fastcore", "fastcore": {"fallback_backend": "scanpy_legacy"}}}
    adata = AnnData(np.ones((4, 3)))
    caps = FastCoreCapabilities(
        omicverse_available=False,
        torch_available=False,
        cuda_available=False,
        rapids_available=False,
        anndataoom_available=False,
        rust_backend_available=False,
        reasons=["omicverse_unavailable"],
    )
    plan = plan_fastcore_backend(cfg, adata=adata, capabilities=caps)
    assert plan == FastCorePlan(
        selected_backend="scanpy_legacy",
        fallback_required=True,
        fallback_backend="scanpy_legacy",
        reasons=["omicverse_unavailable"],
        capabilities=caps.to_dict(),
    )


def test_planner_requires_explicit_experimental_adapter_enable():
    cfg = {"core": {"engine": "fastcore", "fastcore": {"fallback_backend": "scanpy_legacy"}}}
    adata = AnnData(np.ones((4, 3)))
    caps = FastCoreCapabilities(
        omicverse_available=True,
        torch_available=False,
        cuda_available=False,
        rapids_available=False,
        anndataoom_available=False,
        rust_backend_available=False,
        reasons=[],
    )
    plan = plan_fastcore_backend(cfg, adata=adata, capabilities=caps)
    assert plan.selected_backend == "scanpy_legacy"
    assert plan.fallback_required is True
    assert "experimental_omicverse_adapters_disabled" in plan.reasons


def test_planner_selects_cpu_when_enabled_and_available():
    cfg = {
        "core": {
            "engine": "fastcore",
            "fastcore": {
                "fallback_backend": "scanpy_legacy",
                "enable_experimental_omicverse_adapters": True,
                "allowed_backends": ["omicverse_cpu"],
            },
        }
    }
    adata = AnnData(np.ones((4, 3)))
    caps = FastCoreCapabilities(
        omicverse_available=True,
        torch_available=False,
        cuda_available=False,
        rapids_available=False,
        anndataoom_available=False,
        rust_backend_available=False,
        reasons=[],
    )
    plan = plan_fastcore_backend(cfg, adata=adata, capabilities=caps)
    assert plan.selected_backend == "omicverse_cpu"
    assert plan.fallback_required is False
