from __future__ import annotations

import json
from pathlib import Path
from time import perf_counter
from typing import Any, Mapping

from fastcore.backend_plan import canonical_backend, plan_fastcore_backend
from fastcore.backends.scanpy_legacy import run_scanpy_legacy_core
from scsp_agent_sop.assays import apply_final_use, assay_profile
from scsp_agent_sop.config import deep_get
from scsp_agent_sop.core import normalize_log1p
from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.measurement import measurement_record
from scsp_agent_sop.storage import ensure_dir, register_file, write_h5ad_atomic, write_json


MASK_AWARE_BACKENDS = {"scanpy_legacy", "fastcore_cpu"}
# These write log1p_norm themselves (scoring and HVG read it; spatial declares
# its own normalisation). The OmicVerse preprocess backends leave X on their
# shiftlog scale, so the shared definition is written here (S007).
LOG1P_NORM_BACKENDS = {"scanpy_legacy", "fastcore_cpu", "spatial_scanpy"}
# What every core backend hands downstream stages (programmes, annotation, CNV, DE).
CORE_OUTPUTS = {
    "layers": ("log1p_norm",),
    "obsm": ("X_pca_biology", "X_pca_identity_prebatch", "X_pca_harmony_identity", "X_umap_identity", "X_umap_biology"),
    "obsp": ("connectivities_identity",),
    "obs": ("cluster_identity",),
    "var": ("highly_variable_biology",),
}


def missing_core_outputs(adata, counts_layer: str) -> list[str]:
    required = {**CORE_OUTPUTS, "layers": (*CORE_OUTPUTS["layers"], counts_layer)}
    return [f"{slot}[{key!r}]" for slot, keys in required.items() for key in keys if key not in getattr(adata, slot)]


def run_core_pipeline(
    adata,
    cfg: Mapping[str, Any],
    run_root: str | Path,
    *,
    input_path: str | Path | None = None,
    output_path: str | Path | None = None,
) -> dict[str, Any]:
    """Run core through FastCore planning with one allowed legacy fallback."""
    run_root = Path(run_root)
    engine = str(deep_get(cfg, "core.engine", "fastcore"))
    start = perf_counter()
    n_input = apply_final_use(adata) if adata is not None else None
    if assay_profile(cfg) != "scrna":
        if adata is None:
            raise ValueError("Spatial core requires an in-memory AnnData object.")
        plan = None
        engine = "spatial"
        backend = "spatial_scanpy"
        fallback_used = False
        fallback_reason = None
    elif engine == "scanpy_legacy":
        if adata is None:
            raise ValueError("scanpy_legacy requires an in-memory AnnData object.")
        plan = None
        backend = "scanpy_legacy"
        fallback_used = False
        fallback_reason = None
    elif engine == "fastcore":
        plan = plan_fastcore_backend(cfg, adata=adata, input_path=input_path)
        backend = canonical_backend(plan.selected_backend)
        fallback_used = bool(plan.fallback_required)
        fallback_reason = ";".join(plan.reasons) if fallback_used and plan.reasons else None
    else:
        raise ValueError(f"unsupported core.engine: {engine}")

    if backend == "spatial_scanpy":
        from fastcore.backends.spatial_panel import run_spatial_core

        result = run_spatial_core(adata, cfg, run_root)
    elif backend == "scanpy_legacy":
        if adata is None:
            raise ValueError("scanpy_legacy requires an in-memory AnnData object.")
        result = run_scanpy_legacy_core(adata, cfg, run_root)
    elif backend == "fastcore_cpu":
        from fastcore.backends.omicverse_cpu import run_omicverse_cpu_core

        result = run_omicverse_cpu_core(adata, cfg, run_root)
    elif backend == "fastcore_mixed":
        from fastcore.backends.omicverse_mixed import run_omicverse_mixed_core

        result = run_omicverse_mixed_core(adata, cfg, run_root)
    elif backend == "fastcore_oom":
        if input_path is None or output_path is None:
            raise ValueError("fastcore_oom requires input_path and output_path.")
        from fastcore.backends.omicverse_rust_oom import run_omicverse_rust_oom_core

        result = run_omicverse_rust_oom_core(input_path, output_path, cfg, run_root)
    else:
        raise ValueError(f"unsupported FastCore backend: {backend}")

    result = dict(result)
    result_adata = result.pop("_adata", None)
    manifest_adata = result_adata if result_adata is not None else adata
    result["backend"] = backend
    result["engine"] = engine
    result["fallback_used"] = fallback_used
    result["fallback_reason"] = fallback_reason
    counts_layer = str(deep_get(cfg, "qc.counts_layer", "counts"))
    target_sum = float(deep_get(cfg, "core.normalize_total_target_sum", 10000))
    size_factor_key = deep_get(cfg, "core.size_factor_key", None)
    if backend not in LOG1P_NORM_BACKENDS:
        step = perf_counter()
        normalize_log1p(manifest_adata, counts_layer=counts_layer, target_sum=target_sum, size_factor_key=size_factor_key)
        result.setdefault("timings", {})["log1p_norm"] = perf_counter() - step
    log1p_norm = None if backend == "spatial_scanpy" else {
        "target_sum": target_sum, "size_factor_key": size_factor_key, "counts_layer": counts_layer,
        "written_by": backend if backend in LOG1P_NORM_BACKENDS else "core_runner",
    }
    missing = missing_core_outputs(manifest_adata, counts_layer)
    if missing:
        raise RuntimeError(f"{backend} did not deliver the core output contract; missing {', '.join(missing)}")
    result["wall_time_s"] = perf_counter() - start
    # Only these backends select HVGs through the mask-aware select_hvg; the
    # OmicVerse preprocess backends cannot take a mask, so say so instead of
    # letting structural zeros pass as measured.
    result["measurement"] = measurement_record(manifest_adata, honoured=backend in MASK_AWARE_BACKENDS)
    batch_correction = result.get("batch_correction")
    selection = manifest_adata.uns.get("leiden_selection", {}).get("cluster_identity")
    clustering = None
    if selection is not None:
        clustering = {key: value for key, value in selection.items() if key != "candidates_json"}
        clustering["candidates"] = json.loads(selection["candidates_json"])
    review = []
    if result["measurement"]["status"] == "not_honoured":
        review.append(f"{backend} cannot apply the declared feature measurement; HVGs may count unmeasured genes as zeros")
    if batch_correction is not None and batch_correction["status"] == "fallback_uncorrected":
        review.append(f"Harmony was requested but did not run ({batch_correction['reason']}); "
                      "X_pca_harmony_identity is the uncorrected PCA")

    manifest = {
        "schema_version": "fastcore_manifest.v1",
        "qc_input_observations": n_input,
        "engine": engine,
        "backend": result["backend"],
        "gpu": result["backend"] == "fastcore_mixed",
        "n_obs": int(result.get("n_obs", getattr(manifest_adata, "n_obs", 0))),
        "n_vars": int(result.get("n_vars", getattr(manifest_adata, "n_vars", 0))),
        "n_hvgs": int(deep_get(cfg, "core.fastcore.omicverse.n_hvgs", deep_get(cfg, "core.n_top_hvg", 3000))),
        "n_pcs": int(deep_get(cfg, "core.fastcore.omicverse.n_pcs", deep_get(cfg, "core.n_pcs", 50))),
        "timings": result.get("timings", {}),
        "quality": result.get("quality", {}),
        "fallback_used": fallback_used,
        "fallback_reason": fallback_reason,
        "plan": plan.to_dict() if plan is not None else None,
        "artifacts": result.get("artifacts", {}),
        "measurement": result["measurement"],
        "log1p_norm": log1p_norm,
        "batch_correction": batch_correction,
        "clustering": clustering,
    }
    fastcore_dir = ensure_dir(run_root / "core" / "fastcore")
    manifest_path = write_json(manifest, fastcore_dir / "fastcore_manifest.json")
    quality_path = write_json(
        result.get("quality", {"accepted": True, "backend": result["backend"]}), fastcore_dir / "core_quality.json"
    )
    if manifest_adata is not None:
        register_file(
            manifest_adata,
            key="fastcore_manifest",
            path=manifest_path,
            category="artifacts",
            schema="fastcore_manifest.v1",
        )
        register_file(
            manifest_adata, key="core_quality", path=quality_path, category="artifacts", schema="core_quality.v1"
        )
    if result_adata is not None and output_path is not None:
        write_h5ad_atomic(result_adata, output_path)
        result["output_h5ad"] = str(output_path)

    log_decision(
        run_root,
        module="core_analysis",
        decision=f"{engine}_complete",
        reason="Explicit assay routing; RNA retains the FastCore/Scanpy policy and spatial uses its declared preprocessing baseline.",
        parameters={
            "selected_backend": result["backend"],
            "fallback_backend": deep_get(
                cfg, "core.fastcore.fallback_backend", deep_get(cfg, "core.fallback_engine", "scanpy_legacy")
            ),
            "fallback_used": fallback_used,
            "batch_correction_method": result.get(
                "batch_correction_method", deep_get(cfg, "core.batch_correction.method", "harmony2")
            ),
            "harmony2_used": result.get("harmony2_used"),
            "batch_keys": result.get("batch_keys", []),
            "batch_key_source": batch_correction.get("key_source") if batch_correction else None,
            "leiden_selection_policy": clustering["policy"] if clustering else None,
        },
        evidence={
            "n_clusters": result.get("n_clusters"),
            "fastcore_manifest": str(manifest_path),
            "core_quality": str(quality_path),
            "measurement": result["measurement"],
            "batch_correction_status": batch_correction["status"] if batch_correction else None,
            "batch_correction_reason": batch_correction.get("reason") if batch_correction else None,
            "chosen_resolution": clustering["chosen_resolution"] if clustering else None,
        },
        fallback_used=fallback_used,
        fallback_reason=fallback_reason,
        human_review_required=bool(review),
        review_reason="; ".join(review) or None,
    )
    return result
