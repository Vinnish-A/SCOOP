#!/usr/bin/env python
"""Prepare pseudobulk counts and metadata for FastDE/R reference DE.

Usage:
  python scripts/08_prepare_pseudobulk.py --config runs/<run_id>/config/run.yaml \
    --celltype-key cell_type_lvl3 --condition-key condition

The script aggregates raw counts by sample x cell type x condition. It writes
one counts matrix and one metadata table per cell type into 07_de/pseudobulk.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from scsp_agent_sop.config import read_yaml, deep_get, resolve_run_root
from scsp_agent_sop.storage import write_table, register_file, init_file_registry, ensure_dir
from scsp_agent_sop.decision_log import log_decision


def _sum_counts(matrix):
    if sparse.issparse(matrix):
        return np.asarray(matrix.sum(axis=0)).ravel()
    return np.asarray(matrix.sum(axis=0)).ravel()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--celltype-key", default=None)
    ap.add_argument("--condition-key", default=None)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    input_path = Path(args.input) if args.input else run_root / "artifacts" / "final_adata.h5ad"
    if not input_path.exists():
        input_path = run_root / "artifacts" / "adata_annotation_evidence.h5ad"
    adata = ad.read_h5ad(input_path)
    init_file_registry(adata, deep_get(cfg, "run.run_id", run_root.name))

    counts_layer = deep_get(cfg, "qc.counts_layer", "counts")
    sample_key = deep_get(cfg, "keys.sample", "sample_id")
    donor_key = deep_get(cfg, "keys.donor", "donor_id")
    condition_key = args.condition_key or deep_get(cfg, "keys.condition", "condition")
    celltype_key = args.celltype_key or deep_get(cfg, "keys.celltype_primary", "cell_type_lvl3")
    min_cells = int(deep_get(cfg, "de.min_cells_per_sample_celltype", 20))
    for k in [sample_key, condition_key, celltype_key]:
        if k not in adata.obs:
            raise KeyError(f"Missing obs key for pseudobulk: {k}")
    X = adata.layers[counts_layer]
    base = run_root / "07_de" / "pseudobulk"
    ensure_dir(base)
    summary_rows = []
    for ct, idx_ct in adata.obs.groupby(celltype_key).indices.items():
        sub_obs = adata.obs.iloc[list(idx_ct)]
        rows = []
        meta = []
        for key_tuple, idx in sub_obs.groupby([sample_key, condition_key]).indices.items():
            sample, condition = key_tuple
            global_idx = adata.obs.index.get_indexer(sub_obs.iloc[list(idx)].index)
            if len(global_idx) < min_cells:
                continue
            rows.append(_sum_counts(X[global_idx, :]))
            d = {
                "pseudobulk_id": f"{ct}|{sample}|{condition}",
                "cell_type": str(ct),
                sample_key: str(sample),
                condition_key: str(condition),
                "n_cells": len(global_idx),
            }
            if donor_key in adata.obs:
                donor_values = adata.obs.iloc[global_idx][donor_key].astype(str).unique()
                d[donor_key] = donor_values[0] if len(donor_values) == 1 else ";".join(sorted(donor_values))
            meta.append(d)
        if not rows:
            continue
        counts = pd.DataFrame(np.vstack(rows), index=[m["pseudobulk_id"] for m in meta], columns=adata.var_names)
        meta_df = pd.DataFrame(meta)
        safe_ct = str(ct).replace("/", "_").replace(" ", "_")
        ct_dir = base / safe_ct
        c_path = write_table(counts.reset_index(names="pseudobulk_id"), ct_dir / "counts.tsv")
        m_path = write_table(meta_df, ct_dir / "metadata.tsv")
        register_file(adata, key=f"pseudobulk_counts_{safe_ct}", path=c_path, schema="pseudobulk_counts.v1")
        register_file(adata, key=f"pseudobulk_metadata_{safe_ct}", path=m_path, schema="pseudobulk_metadata.v1")
        summary_rows.append({"cell_type": ct, "n_pseudobulk_samples": len(meta_df), "total_cells": int(meta_df["n_cells"].sum())})
    summary = pd.DataFrame(summary_rows)
    sum_path = write_table(summary, base / "pseudobulk_summary.tsv")
    register_file(adata, key="pseudobulk_summary", path=sum_path, schema="pseudobulk_summary.v1")
    adata.write_h5ad(run_root / "artifacts" / "adata_pseudobulk_indexed.h5ad")
    log_decision(run_root, module="de", decision="pseudobulk_prepared", reason="Raw counts aggregated by sample x condition x cell type for FastDE condition DE and R reference validation.", parameters={"celltype_key": celltype_key, "condition_key": condition_key, "min_cells": min_cells}, evidence={"n_celltypes": len(summary_rows)})


if __name__ == "__main__":
    main()
