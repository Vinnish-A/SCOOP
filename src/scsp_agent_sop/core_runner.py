from __future__ import annotations

from pathlib import Path
from time import perf_counter
from typing import Any, Mapping

from fastcore.backend_plan import plan_fastcore_backend
from fastcore.backends.scanpy_legacy import run_scanpy_legacy_core
from scsp_agent_sop.config import deep_get
from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.storage import ensure_dir, register_file, write_json


def run_core_pipeline(
    adata,
    cfg: Mapping[str, Any],
    run_root: str | Path,
    *,
    input_path: str | Path | None = None,
    output_path: str | Path | None = None,
) -> dict[str, Any]:
    """Run 02_core through FastCore planning with one allowed legacy fallback."""
    run_root = Path(run_root)
    engine = str(deep_get(cfg, "core.engine", "fastcore"))
    start = perf_counter()
    if engine == "scanpy_legacy":
        plan = None
        backend = "scanpy_legacy"
        fallback_used = False
        fallback_reason = None
    elif engine == "fastcore":
        plan = plan_fastcore_backend(cfg, adata=adata, input_path=input_path)
        backend = plan.selected_backend
        fallback_used = bool(plan.fallback_required)
        fallback_reason = ";".join(plan.reasons) if fallback_used and plan.reasons else None
    else:
        raise ValueError(f"unsupported core.engine: {engine}")

    if backend == "scanpy_legacy":
        result = run_scanpy_legacy_core(adata, cfg, run_root)
    elif backend == "omicverse_cpu":
        from fastcore.backends.omicverse_cpu import run_omicverse_cpu_core

        result = run_omicverse_cpu_core(adata, cfg, run_root)
    elif backend == "omicverse_cpu_gpu_mixed":
        from fastcore.backends.omicverse_mixed import run_omicverse_mixed_core

        result = run_omicverse_mixed_core(adata, cfg, run_root)
    elif backend == "omicverse_gpu_rapids":
        from fastcore.backends.omicverse_gpu import run_omicverse_gpu_core

        result = run_omicverse_gpu_core(adata, cfg, run_root)
    elif backend == "omicverse_rust_oom":
        if input_path is None or output_path is None:
            raise ValueError("omicverse_rust_oom requires input_path and output_path.")
        from fastcore.backends.omicverse_rust_oom import run_omicverse_rust_oom_core

        result = run_omicverse_rust_oom_core(input_path, output_path, cfg, run_root)
    else:
        raise ValueError(f"unsupported FastCore backend: {backend}")

    result = dict(result)
    result.setdefault("backend", backend)
    result["engine"] = engine
    result["fallback_used"] = fallback_used
    result["fallback_reason"] = fallback_reason
    result["wall_time_s"] = perf_counter() - start

    manifest = {
        "schema_version": "fastcore_manifest.v1",
        "engine": engine,
        "backend": result["backend"],
        "gpu": result["backend"] in {"omicverse_cpu_gpu_mixed", "omicverse_gpu_rapids"},
        "rapids": result["backend"] == "omicverse_gpu_rapids",
        "n_obs": int(result.get("n_obs", adata.n_obs)),
        "n_vars": int(result.get("n_vars", adata.n_vars)),
        "n_hvgs": int(deep_get(cfg, "core.fastcore.omicverse.n_hvgs", deep_get(cfg, "core.n_top_hvg", 3000))),
        "n_pcs": int(deep_get(cfg, "core.fastcore.omicverse.n_pcs", deep_get(cfg, "core.n_pcs", 50))),
        "timings": result.get("timings", {}),
        "quality": result.get("quality", {}),
        "fallback_used": fallback_used,
        "fallback_reason": fallback_reason,
        "plan": plan.to_dict() if plan is not None else None,
        "artifacts": result.get("artifacts", {}),
    }
    fastcore_dir = ensure_dir(run_root / "02_core" / "fastcore")
    manifest_path = write_json(manifest, fastcore_dir / "fastcore_manifest.json")
    quality_path = write_json(result.get("quality", {"accepted": True, "backend": result["backend"]}), fastcore_dir / "core_quality.json")
    register_file(adata, key="fastcore_manifest", path=manifest_path, category="artifacts", schema="fastcore_manifest.v1")
    register_file(adata, key="core_quality", path=quality_path, category="artifacts", schema="core_quality.v1")

    log_decision(
        run_root,
        module="core_analysis",
        decision="fastcore_complete" if engine == "fastcore" else "scanpy_legacy_complete",
        reason="02_core executed through FastCore planner; Scanpy legacy is the only non-FastCore fallback.",
        parameters={
            "selected_backend": result["backend"],
            "fallback_backend": deep_get(cfg, "core.fastcore.fallback_backend", deep_get(cfg, "core.fallback_engine", "scanpy_legacy")),
            "fallback_used": fallback_used,
            "harmony_pytorch_used": result.get("harmony_pytorch_used"),
            "batch_keys": result.get("batch_keys", []),
        },
        evidence={
            "n_clusters": result.get("n_clusters"),
            "fastcore_manifest": str(manifest_path),
            "core_quality": str(quality_path),
        },
        fallback_used=fallback_used,
        fallback_reason=fallback_reason,
    )
    return result
