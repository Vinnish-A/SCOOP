#!/usr/bin/env python
"""Run FastCNMF-preferred programme discovery.

Usage:
  scoop programs --config runs/<run_id>/config/run.yaml

The native stage reports its actual sklearn multistart/consensus implementation.
Legacy method names remain accepted as compatibility aliases. The standalone
FastCNMF executable is a separate implementation, not silently claimed here.
With programs.batch_correction: harmony_moe, FastCNMF preprocessing corrects
seurat_v3 HVGs for the declared run.obs_keys.batch before the same consensus NMF
(S010); a guessed keys.batch_candidates column is refused. A declared feature
measurement is recorded as not_honoured: NMF reads unmeasured genes as zeros.

"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import anndata as ad
import pandas as pd

from scsp_agent_sop.config import deep_get, read_yaml, resolve_obs_keys, resolve_run_root
from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.measurement import measurement_record
from scsp_agent_sop.programs import (
    FASTCNMF_COMPATIBLE_METHODS,
    FASTCNMF_DEFAULT_MAX_ITER,
    FASTCNMF_DEFAULT_SEEDS,
    run_fast_consensus_nmf,
    run_harmony_consensus_nmf,
    run_omicverse_cnmf_validation,
)
from scsp_agent_sop.storage import write_h5ad_atomic, ensure_dir, init_file_registry, register_file, write_table


def main() -> dict:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--output", default=None)
    ap.add_argument("--validate-with-omicverse-cnmf", action="store_true")
    args = ap.parse_args()
    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    input_path = Path(args.input) if args.input else run_root / "artifacts" / "adata_core.h5ad"
    output = Path(args.output) if args.output else run_root / "artifacts" / "adata_programs.h5ad"
    adata = ad.read_h5ad(input_path)
    init_file_registry(adata, deep_get(cfg, "run.run_id", run_root.name))

    method = str(deep_get(cfg, "programs.method", "sklearn_multistart_nmf")).lower()
    if method not in FASTCNMF_COMPATIBLE_METHODS:
        supported = ", ".join(sorted(FASTCNMF_COMPATIBLE_METHODS))
        raise ValueError(f"unsupported programs.method={method!r}; supported methods: {supported}")

    seeds = list(deep_get(cfg, "programs.seeds", list(FASTCNMF_DEFAULT_SEEDS)))
    max_iter = int(deep_get(cfg, "programs.max_iter", FASTCNMF_DEFAULT_MAX_ITER))
    batch_correction = str(deep_get(cfg, "programs.batch_correction", "none"))
    if batch_correction not in {"none", "harmony_moe"}:
        raise ValueError(f"unsupported programs.batch_correction={batch_correction!r}; use none or harmony_moe")
    keys = resolve_obs_keys(cfg, adata.obs, required=("batch",) if batch_correction == "harmony_moe" else ())
    if batch_correction == "harmony_moe" and keys.sources["batch"] != "run.obs_keys.batch":
        # Expression-level correction regresses out whatever it is given; a
        # sample or guessed candidate column is not a declared technical batch.
        raise ValueError(f"programs.batch_correction: harmony_moe needs a declared run.obs_keys.batch; "
                         f"{keys.batch!r} came from {keys.sources['batch']}, a candidate guess")
    # Neither NMF route restricts genes per sample, so a declared measurement is not honoured.
    measurement = measurement_record(adata, honoured=False)
    nmf_kwargs = dict(
        k_grid=deep_get(cfg, "programs.k_grid", [5, 8, 10, 12, 15]),
        seeds=seeds,
        max_iter=max_iter,
        workers=int(deep_get(cfg, "programs.workers", 4)),
        stability_threshold=deep_get(cfg, "programs.stability_threshold", 0.70),
    )
    preprocess = None
    if batch_correction == "harmony_moe":
        summary, weights, usage, preprocess = run_harmony_consensus_nmf(
            adata,
            input_h5ad=input_path,
            batch_key=keys.batch,
            output_prefix=run_root / "programs" / "fastcnmf_preprocess" / "programs",
            n_top_genes=int(deep_get(cfg, "programs.n_top_hvg", 3000)),
            seed=int(deep_get(cfg, "run.random_seed", 0)),
            **nmf_kwargs,
        )
    else:
        summary, weights, usage = run_fast_consensus_nmf(
            adata,
            layer=deep_get(cfg, "programs.layer", "log1p_norm"),
            hvg_key=deep_get(cfg, "programs.hvg_key", "highly_variable_biology"),
            **nmf_kwargs,
        )
    batch_handling = {"method": batch_correction, **keys.record()["batch"]}
    summary_path = write_table(summary, run_root / "programs" / "tables" / "nmf_k_sweep.tsv")
    weights_path = write_table(weights, run_root / "programs" / "tables" / "nmf_gene_weights.parquet")
    usage_path = write_table(usage, run_root / "programs" / "tables" / "nmf_usage.parquet")
    stability_path = write_table(pd.DataFrame(adata.uns["nmf_component_evidence"]), run_root / "programs/tables/nmf_component_stability.tsv")
    register_file(adata, key="nmf_k_sweep", path=summary_path, schema="nmf_k_sweep.v1")
    register_file(adata, key="nmf_gene_weights", path=weights_path, schema="nmf_gene_weights.v1")
    register_file(adata, key="nmf_usage", path=usage_path, schema="nmf_usage.v1")
    register_file(adata, key="nmf_component_stability", path=stability_path, schema="nmf_component_stability.v2")

    selected_k = int(summary.loc[summary["selected"], "k"].iloc[0])
    min_stability = float(summary.loc[summary.selected, "min_component_stability"].iloc[0])
    convergence_fraction = float(summary.loc[summary.selected, "converged_fraction"].iloc[0])
    threshold = float(deep_get(cfg, "programs.stability_threshold", 0.70))
    review_required = (not math.isfinite(min_stability) or min_stability < threshold or convergence_fraction < 1.0
                       or not adata.uns["nmf_method"]["usage_refit_converged"])
    fallback_used = False
    fallback_reason = None
    if args.validate_with_omicverse_cnmf or deep_get(cfg, "programs.validate_with_omicverse_cnmf", False):
        try:
            run_omicverse_cnmf_validation(
                adata,
                components=deep_get(cfg, "programs.k_grid", [5, 8, 10, 12, 15]),
                output_dir=str(run_root / "programs" / "omicverse_cnmf_validation"),
                n_iter=100,
                use_gpu=True,
            )
            fallback_used = True
            fallback_reason = "OmicVerse cNMF validation explicitly requested"
        except Exception as exc:
            fallback_used = True
            fallback_reason = f"OmicVerse cNMF validation attempted but failed: {type(exc).__name__}: {exc}"

    review_reasons = []
    if review_required:
        review_reasons.append("NMF stability is low/unassessed or some selected-K replicates did not converge; "
                              "inspect programme tables.")
    if measurement["status"] == "not_honoured":
        review_reasons.append("Programme discovery cannot apply the declared feature measurement; "
                              "NMF reads unmeasured genes as zeros.")
    ensure_dir(output.parent)
    write_h5ad_atomic(adata, output)
    log_decision(
        run_root,
        module="programs",
        decision="sklearn_multistart_nmf_complete",
        reason="Random-start sklearn NMF, component consensus and fixed-component usage refit; convergence reported.",
        parameters={
            "requested_method": method,
            "actual_method": "sklearn_multistart_nmf",
            "selected_k": selected_k,
            "k_grid": deep_get(cfg, "programs.k_grid", []),
            "n_seeds": len(seeds),
            "max_iter": max_iter,
            "batch_handling": batch_handling,
        },
        evidence={"min_stability": min_stability, "selected_converged_fraction": convergence_fraction,
                  "measurement": measurement},
        fallback_used=fallback_used,
        fallback_reason=fallback_reason,
        human_review_required=bool(review_reasons),
        review_reason="; ".join(review_reasons) or None,
    )
    limitations = ["Discovery feature coverage differs from whole-input projection; programme usage does not assign identity."]
    if preprocess is not None:
        limitations.append(f"Programmes were discovered on {preprocess['n_vars']} Harmony-corrected seurat_v3 HVGs "
                           f"(batch {keys.batch!r}); the HVGs were chosen batch-blind and the matrix is dense (S012).")
    elif keys.sources["batch"] == "run.obs_keys.batch":
        limitations.append(f"Programme discovery did not model the declared batch {keys.batch!r}; "
                           "set programs.batch_correction: harmony_moe to correct for it.")
    if measurement["status"] == "not_honoured":
        limitations.append(f"The declared feature measurement ({measurement['path']}) was not applied; genes a sample "
                           "did not measure enter NMF as zeros.")
    preprocess_outputs = [preprocess["corrected_h5ad"], preprocess["hvg_txt"],
                          str(Path(preprocess["output_prefix"]).parent / "fastcnmf_preprocess_manifest.json")] if preprocess else []
    return {"inputs": [str(input_path.resolve())],
            "outputs": [str(output.resolve()), str(summary_path.resolve()), str(weights_path.resolve()), str(usage_path.resolve()), str(stability_path.resolve()),
                        *preprocess_outputs],
            "execution": {"status": "completed"},
            "validation": {"status": "limited" if review_required else "passed", "scope": "cross_seed_component_stability_and_fit_termination",
                           "checks": [{"name": "selected_components_cross_seed_stability_and_convergence", "status": "limited" if review_required else "passed"}]},
            "interpretation": {"status": "not_assessed"},
            "coverage": [{"scope": "programme_usage", "unit": "observation_views", "total": adata.n_obs, "processed": adata.n_obs}],
            "limitations": limitations,
            "programme_method": dict(adata.uns["nmf_method"]),
            "batch_handling": batch_handling,
            "measurement": measurement}


if __name__ == "__main__":
    main()
