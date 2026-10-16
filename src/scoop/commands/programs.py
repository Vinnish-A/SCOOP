#!/usr/bin/env python
"""Run FastCNMF-preferred programme discovery.

Usage:
  scoop programs --config runs/<run_id>/config/run.yaml

The native stage reports its actual sklearn multistart/consensus implementation.
Legacy method names remain accepted as compatibility aliases. The standalone
FastCNMF executable is a separate implementation, not silently claimed here.

"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import anndata as ad

from scsp_agent_sop.config import deep_get, read_yaml, resolve_run_root
from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.programs import (
    FASTCNMF_COMPATIBLE_METHODS,
    FASTCNMF_DEFAULT_MAX_ITER,
    FASTCNMF_DEFAULT_SEEDS,
    run_fast_consensus_nmf,
    run_omicverse_cnmf_validation,
)
from scsp_agent_sop.storage import write_h5ad_atomic, ensure_dir, init_file_registry, register_file, write_table


def main() -> None:
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
    summary, weights, usage = run_fast_consensus_nmf(
        adata,
        layer=deep_get(cfg, "programs.layer", "log1p_norm"),
        k_grid=deep_get(cfg, "programs.k_grid", [5, 8, 10, 12, 15]),
        seeds=seeds,
        max_iter=max_iter,
        workers=int(deep_get(cfg, "programs.workers", 4)),
        stability_threshold=deep_get(cfg, "programs.stability_threshold", 0.70),
    )
    summary_path = write_table(summary, run_root / "programs" / "tables" / "nmf_k_sweep.tsv")
    weights_path = write_table(weights, run_root / "programs" / "tables" / "nmf_gene_weights.parquet")
    usage_path = write_table(usage, run_root / "programs" / "tables" / "nmf_usage.parquet")
    register_file(adata, key="nmf_k_sweep", path=summary_path, schema="nmf_k_sweep.v1")
    register_file(adata, key="nmf_gene_weights", path=weights_path, schema="nmf_gene_weights.v1")
    register_file(adata, key="nmf_usage", path=usage_path, schema="nmf_usage.v1")

    selected_k = int(summary.loc[summary["selected"], "k"].iloc[0])
    min_stability = float(summary["median_stability"].min())
    convergence_fraction = float(summary.loc[summary.selected, "converged_fraction"].iloc[0])
    review_required = not math.isfinite(min_stability) or min_stability < 0.70 or convergence_fraction < 1.0
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
        },
        evidence={"min_stability": min_stability, "selected_converged_fraction": convergence_fraction},
        fallback_used=fallback_used,
        fallback_reason=fallback_reason,
        human_review_required=review_required,
        review_reason=(
            "NMF stability is low/unassessed or some selected-K replicates did not converge; inspect programme tables."
        )
        if review_required
        else None,
    )


if __name__ == "__main__":
    main()
