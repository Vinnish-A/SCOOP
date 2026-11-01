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


def test_planner_reads_h5ad_size_and_does_not_treat_a_path_as_backed(tmp_path):
    path = tmp_path / "input.h5ad"
    AnnData(np.ones((4, 3))).write_h5ad(path)
    caps = FastCoreCapabilities(
        omicverse_available=True,
        torch_available=False,
        cuda_available=False,
        anndataoom_available=True,
        rust_backend_available=True,
        vendored_omicverse_available=True,
        reasons=[],
    )
    cfg = {"core": {"fastcore": {}}}
    assert plan_fastcore_backend(cfg, input_path=path, capabilities=caps).selected_backend == "fastcore_cpu"
    declared = {"core": {"fastcore": {"auto": {"input_backed": True}}}}
    assert plan_fastcore_backend(declared, input_path=path, capabilities=caps).selected_backend == "fastcore_oom"
    large = {"core": {"fastcore": {"auto": {"rust_oom_min_cells": 4}}}}
    assert plan_fastcore_backend(large, input_path=path, capabilities=caps).selected_backend == "fastcore_oom"


def test_planner_sizes_an_h5ad_from_obs_and_the_counts_layer_without_x(tmp_path):
    import pandas as pd
    from scipy import sparse

    path = tmp_path / "layers_only.h5ad"
    adata = AnnData(obs=pd.DataFrame(index=[f"c{i}" for i in range(5)]), var=pd.DataFrame(index=list("abc")))
    adata.layers["counts"] = sparse.csr_matrix(np.ones((5, 3), np.float32))  # nnz 15; X is absent
    adata.write_h5ad(path)
    caps = FastCoreCapabilities(
        omicverse_available=True,
        torch_available=False,
        cuda_available=False,
        anndataoom_available=True,
        rust_backend_available=True,
        vendored_omicverse_available=True,
        reasons=[],
    )
    auto = {"rust_oom_min_cells": 10**6, "rust_oom_min_nnz": 15}
    assert plan_fastcore_backend({"core": {"fastcore": {"auto": auto}}}, input_path=path,
                                 capabilities=caps).selected_backend == "fastcore_oom"
    auto["rust_oom_min_nnz"] = 16
    assert plan_fastcore_backend({"core": {"fastcore": {"auto": auto}}}, input_path=path,
                                 capabilities=caps).selected_backend == "fastcore_cpu"
    assert plan_fastcore_backend({"core": {"fastcore": {"auto": {"rust_oom_min_cells": 5}}}}, input_path=path,
                                 capabilities=caps).selected_backend == "fastcore_oom"
