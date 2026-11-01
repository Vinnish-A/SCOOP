#!/usr/bin/env python
"""Test a v2 reference package's identifiability on holdout pseudo-spots (F12).

Usage:
  scoop spatial identifiability --reference ref/ --output ref/identifiability [--preset standard --threads 4]

The fit uses the managed RCTD adapter with the package's held-out genes excluded. Pass the
written component_identifiability.parquet as spatial.rctd.identifiability_table so that each
fit's component summary carries the verdicts.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from scsp_agent_sop.spatial_identifiability import run_identifiability


def main() -> dict:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True, help="scoop.reference-profiles.v2 package directory")
    ap.add_argument("--output", required=True)
    ap.add_argument("--preset", default="standard", choices=["fast", "standard"])
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--q-matrices", help="RCTD likelihood lookup; default: the backend's cached copy")
    ap.add_argument("--run-purpose", default="pilot", choices=["production", "sensitivity", "pilot"])
    ap.add_argument("--spots-per-component", type=int, default=40)
    ap.add_argument("--dominant-cells", type=int, default=3)
    ap.add_argument("--other-cells", type=int, default=3)
    ap.add_argument("--depth-median", type=float, default=4000)
    ap.add_argument("--depth-sdlog", type=float, default=0.6)
    ap.add_argument("--seed", type=int, default=20261031)
    ap.add_argument("--config", help="run config; only places the stage receipt under its run root")
    args = ap.parse_args()
    options = {"threads": args.threads, **({"q_matrices": args.q_matrices} if args.q_matrices else {})}
    summary = run_identifiability(args.reference, args.output, rctd_options=options, run_purpose=args.run_purpose,
                                  preset=args.preset, spots_per_component=args.spots_per_component,
                                  dominant_cells=args.dominant_cells, other_cells=args.other_cells,
                                  depth_median=args.depth_median, depth_sdlog=args.depth_sdlog, seed=args.seed)
    out = Path(args.output).resolve()
    print(out / "identifiability_summary.json")
    return {"inputs": [str(Path(args.reference).resolve() / "reference_manifest.json")],
            "outputs": [str(out / name) for name in ["identifiability_summary.json", "component_identifiability.parquet",
                                                     "pair_identifiability.parquet"]],
            "execution": {"status": "completed"},
            "validation": {"status": "limited", "scope": "holdout_pseudospots_without_platform_effect",
                           "checks": [{"name": "pseudospots_fitted", "status": "passed" if summary["n_fitted"] else "failed"}]},
            "interpretation": {"status": "not_applicable"},
            "coverage": [{"scope": "holdout_pseudospots", "unit": "pseudo_spots", "total": summary["n_pseudospots"],
                          "processed": summary["n_fitted"]}],
            "limitations": [summary["caveat"]]}


if __name__ == "__main__":
    main()
