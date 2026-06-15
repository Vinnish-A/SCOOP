#!/usr/bin/env python
"""Run FastCNMF-preferred programme discovery.

Usage:
  python scripts/03_fast_consensus_nmf.py --config runs/<run_id>/config/run.yaml

The default method is the SCOOP FastCNMF profile: exact coordinate-descent NMF,
20 replicate seeds and max_iter=50. Legacy sklearn consensus NMF method names
remain accepted for existing configs. OmicVerse cNMF is an explicit validation
fallback for unstable or claim-critical programmes.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import anndata as ad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from scsp_agent_sop.config import read_yaml, deep_get, resolve_run_root
from scsp_agent_sop.programs import (
    FASTCNMF_COMPATIBLE_METHODS,
    FASTCNMF_DEFAULT_MAX_ITER,
    FASTCNMF_DEFAULT_SEEDS,
    run_fast_consensus_nmf,
    run_omicverse_cnmf_validation,
)
from scsp_agent_sop.storage import write_table, register_file, init_file_registry, ensure_dir
from scsp_agent_sop.decision_log import log_decision


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

    method = str(deep_get(cfg, "programs.method", "fastcnmf")).lower()
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
        stability_threshold=deep_get(cfg, "programs.stability_threshold", 0.70),
    )
    summary_path = write_table(summary, run_root / "03_programs" / "tables" / "nmf_k_sweep.tsv")
    weights_path = write_table(weights, run_root / "03_programs" / "tables" / "nmf_gene_weights.parquet")
    usage_path = write_table(usage, run_root / "03_programs" / "tables" / "nmf_usage.parquet")
    register_file(adata, key="nmf_k_sweep", path=summary_path, schema="nmf_k_sweep.v1")
    register_file(adata, key="nmf_gene_weights", path=weights_path, schema="nmf_gene_weights.v1")
    register_file(adata, key="nmf_usage", path=usage_path, schema="nmf_usage.v1")

    selected_k = int(summary.loc[summary["selected"], "k"].iloc[0])
    min_stability = float(summary["median_stability"].min())
    fallback_used = False
    fallback_reason = None
    if args.validate_with_omicverse_cnmf or min_stability < deep_get(cfg, "programs.stability_threshold", 0.70):
        try:
            run_omicverse_cnmf_validation(
                adata,
                components=deep_get(cfg, "programs.k_grid", [5, 8, 10, 12, 15]),
                output_dir=str(run_root / "03_programs" / "omicverse_cnmf_validation"),
                n_iter=100,
                use_gpu=True,
            )
            fallback_used = True
            fallback_reason = "OmicVerse cNMF validation requested or fast NMF stability was below threshold"
        except Exception as exc:
            fallback_used = True
            fallback_reason = f"OmicVerse cNMF validation attempted but failed: {type(exc).__name__}: {exc}"

    ensure_dir(output.parent)
    adata.write_h5ad(output)
    log_decision(
        run_root,
        module="programs",
        decision=f"{method}_complete",
        reason="FastCNMF-preferred programme discovery executed and gene weights externalized.",
        parameters={
            "method": method,
            "selected_k": selected_k,
            "k_grid": deep_get(cfg, "programs.k_grid", []),
            "n_seeds": len(seeds),
            "max_iter": max_iter,
        },
        evidence={"min_stability": min_stability},
        fallback_used=fallback_used,
        fallback_reason=fallback_reason,
        human_review_required=min_stability < 0.70,
        review_reason="At least one K has low consensus stability; inspect programme tables and consider lineage-specific NMF." if min_stability < 0.70 else None,
    )


if __name__ == "__main__":
    main()
