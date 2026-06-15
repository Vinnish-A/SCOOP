#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import resource
import sys
from pathlib import Path
from time import perf_counter

import anndata as ad

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from scsp_agent_sop.config import read_yaml, resolve_run_root
from scsp_agent_sop.core_runner import run_core_pipeline
from scsp_agent_sop.storage import ensure_dir, init_file_registry


def _peak_rss_mb() -> float:
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return float(usage / 1024)


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark one planned FastCore backend on an input H5AD.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-h5ad", default=None)
    parser.add_argument("--output-json", default=None)
    args = parser.parse_args()

    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    input_path = Path(args.input)
    output_h5ad = Path(args.output_h5ad) if args.output_h5ad else run_root / "artifacts" / "adata_core_benchmark.h5ad"
    output_json = Path(args.output_json) if args.output_json else run_root / "02_core" / "fastcore" / "benchmark_core_backends.json"

    adata = ad.read_h5ad(input_path)
    init_file_registry(adata, str(cfg.get("run", {}).get("run_id", run_root.name)))
    start = perf_counter()
    result = run_core_pipeline(adata, cfg, run_root, input_path=input_path, output_path=output_h5ad)
    wall_time_s = perf_counter() - start
    ensure_dir(output_h5ad.parent)
    adata.write_h5ad(output_h5ad)

    record = {
        "backend": result.get("backend"),
        "fallback_used": result.get("fallback_used"),
        "fallback_reason": result.get("fallback_reason"),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "wall_time_s": wall_time_s,
        "peak_rss_mb": _peak_rss_mb(),
        "timings": result.get("timings", {}),
        "quality": result.get("quality", {}),
        "output_h5ad": str(output_h5ad),
    }
    ensure_dir(output_json.parent)
    output_json.write_text(json.dumps(record, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(record, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
