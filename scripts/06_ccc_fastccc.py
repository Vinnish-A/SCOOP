#!/usr/bin/env python
"""Run FastCCC screening and optional complex-aware validation.

Usage:
  python scripts/06_ccc_fastccc.py --config runs/<run_id>/config/run.yaml
  python scripts/06_ccc_fastccc.py --config runs/<run_id>/config/run.yaml --dry-run

FastCCC output is externalized. CellPhoneDB/LIANA validation is optional and
only meant for complex-sensitive or claim-critical interactions.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import anndata as ad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from scsp_agent_sop.config import read_yaml, deep_get, resolve_run_root, project_root_from_run_root
from scsp_agent_sop.ccc import validate_lr_resource, choose_ccc_groupby, build_fastccc_command, run_fastccc_command, run_cellphonedb_validation_omicverse, run_liana_validation_omicverse
from scsp_agent_sop.storage import write_table, write_json, register_file, init_file_registry, ensure_dir
from scsp_agent_sop.decision_log import log_decision


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--validate-cellphonedb", action="store_true")
    ap.add_argument("--validate-liana", action="store_true")
    args = ap.parse_args()
    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    project_root = project_root_from_run_root(run_root)
    if not deep_get(cfg, "ccc.enabled", False):
        print("CCC module disabled in config. Set ccc.enabled=true to run.")
        return
    input_h5ad = Path(args.input) if args.input else run_root / "artifacts" / "adata_annotation_evidence.h5ad"
    adata = ad.read_h5ad(input_h5ad)
    init_file_registry(adata, deep_get(cfg, "run.run_id", run_root.name))

    lr_resource = project_root / deep_get(cfg, "ccc.lr_resource")
    lr = validate_lr_resource(lr_resource)
    complex_path = write_table(lr[lr["complex_sensitive"]].copy(), run_root / "06_ccc" / "fastccc" / "tables" / "complex_sensitive_lris.parquet")
    register_file(adata, key="complex_sensitive_lris", path=complex_path, schema="lr_resource_complex_sensitive.v1")

    groupby = choose_ccc_groupby(
        adata,
        primary=deep_get(cfg, "ccc.groupby", "cell_type_lvl3"),
        fallback=deep_get(cfg, "ccc.fallback_groupby", "cell_type_lvl2"),
        min_cells_per_group=deep_get(cfg, "ccc.min_cells_per_group", 20),
    )
    output_dir = run_root / "06_ccc" / "fastccc"
    ensure_dir(output_dir / "tables")
    cmd = build_fastccc_command(
        deep_get(cfg, "ccc.fastccc_command_template"),
        input_h5ad=input_h5ad,
        groupby=groupby,
        lr_resource=lr_resource,
        output_dir=output_dir,
    )
    write_json({"command": cmd, "groupby": groupby}, run_root / "06_ccc" / "fastccc" / "fastccc_command.json")
    if args.dry_run:
        print(" ".join(cmd))
        return
    completed = run_fastccc_command(cmd, cwd=project_root)

    fallback_used = False
    fallback_reason = None
    if args.validate_cellphonedb or deep_get(cfg, "ccc.complex_validation.enabled", False):
        try:
            val = run_cellphonedb_validation_omicverse(
                adata,
                cpdb_file_path=str(project_root / deep_get(cfg, "ccc.complex_validation.cpdb_file_path")),
                celltype_key=groupby,
                output_dir=run_root / "06_ccc" / "cellphonedb" / "tables",
                iterations=deep_get(cfg, "ccc.complex_validation.iterations", 1000),
                threshold=deep_get(cfg, "ccc.complex_validation.threshold", 0.1),
                pvalue=deep_get(cfg, "ccc.complex_validation.pvalue", 0.05),
                threads=deep_get(cfg, "ccc.complex_validation.threads", 10),
            )
            write_json(val, run_root / "06_ccc" / "cellphonedb" / "cellphonedb_export_summary.json")
        except Exception as exc:
            fallback_used = True
            fallback_reason = f"CellPhoneDB validation failed: {type(exc).__name__}: {exc}"
    if args.validate_liana or deep_get(cfg, "ccc.liana_validation.enabled", False):
        try:
            liana = run_liana_validation_omicverse(adata, groupby=groupby, method=deep_get(cfg, "ccc.liana_validation.method", "rank_aggregate"))
            p = write_table(liana, run_root / "06_ccc" / "liana" / "tables" / "liana_rank_aggregate.parquet")
            register_file(adata, key="liana_rank_aggregate", path=p, schema="liana_rank_aggregate.v1")
        except Exception as exc:
            fallback_used = True
            fallback_reason = (fallback_reason or "") + f"; LIANA validation failed: {type(exc).__name__}: {exc}"

    adata.uns.setdefault("ccc", {})
    adata.uns["ccc"].update({"primary_method": "FastCCC", "groupby": groupby, "lr_resource": str(lr_resource), "complex_sensitive_count": int(lr["complex_sensitive"].sum())})
    output_h5ad = run_root / "artifacts" / "adata_ccc_indexed.h5ad"
    adata.write_h5ad(output_h5ad)
    log_decision(
        run_root,
        module="ccc",
        decision="fastccc_executed",
        reason="FastCCC used as primary screening; complex-sensitive LRIs externalized for validation.",
        parameters={"groupby": groupby, "command": cmd},
        evidence={"complex_sensitive_lris": int(lr['complex_sensitive'].sum()), "stdout_tail": completed.stdout[-1000:]},
        fallback_used=fallback_used,
        fallback_reason=fallback_reason,
    )


if __name__ == "__main__":
    main()
