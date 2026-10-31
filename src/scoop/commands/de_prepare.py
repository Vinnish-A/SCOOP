#!/usr/bin/env python
"""Prepare pseudobulk counts and metadata for FastDE/R reference DE.

Usage:
  scoop de prepare --config runs/<run_id>/config/run.yaml \
    --celltype-key cell_type_lvl3 --condition-key condition

The script aggregates raw counts by sample x cell type x condition. It writes
one counts matrix and one metadata table per cell type into differential-expression/pseudobulk.
With committed annotation, only accepted observations are aggregated (S017). The
metadata carries donor, batch, study and keys.batch_candidates columns for the
DE design.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from scsp_agent_sop.config import deep_get, read_yaml, resolve_obs_keys, resolve_run_root
from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.storage import ensure_dir, init_file_registry, register_file, write_table


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
    keys = resolve_obs_keys(cfg, adata.obs, required=("sample",))
    sample_key = keys.sample
    condition_key = args.condition_key or deep_get(cfg, "keys.condition", "condition")
    celltype_key = args.celltype_key or deep_get(cfg, "keys.celltype_primary", "cell_type_lvl3")
    min_cells = int(deep_get(cfg, "de.min_cells_per_sample_celltype", 20))
    for k in [condition_key, celltype_key]:
        if k not in adata.obs:
            raise KeyError(f"Missing obs key for pseudobulk: {k}")
    # Unreviewed, unresolved and technically excluded cells are not a cell type's sample (S017).
    status_filtered = "annotation_status" in adata.obs
    obs = adata.obs[adata.obs.annotation_status.astype(str).isin(["accepted", "labeled"])] if status_filtered else adata.obs
    design_keys = [k for k in dict.fromkeys([keys.donor, keys.batch, keys.study, *deep_get(cfg, "keys.batch_candidates", [])])
                   if k and k in obs and k not in {sample_key, condition_key}]
    X = adata.layers[counts_layer]
    base = run_root / "differential-expression" / "pseudobulk"
    ensure_dir(base)
    summary_rows = []
    for ct, idx_ct in obs.groupby(celltype_key).indices.items():
        sub_obs = obs.iloc[list(idx_ct)]
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
            for key in design_keys:
                values = adata.obs.iloc[global_idx][key].astype(str).unique()
                d[key] = values[0] if len(values) == 1 else ";".join(sorted(values))
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
        summary_rows.append(
            {"cell_type": ct, "n_pseudobulk_samples": len(meta_df), "total_cells": int(meta_df["n_cells"].sum())}
        )
    summary = pd.DataFrame(summary_rows)
    sum_path = write_table(summary, base / "pseudobulk_summary.tsv")
    register_file(adata, key="pseudobulk_summary", path=sum_path, schema="pseudobulk_summary.v1")
    adata.write_h5ad(run_root / "artifacts" / "adata_pseudobulk_indexed.h5ad")
    log_decision(
        run_root,
        module="de",
        decision="pseudobulk_prepared",
        reason=(
            "Raw counts aggregated by sample x condition x cell type for FastDE "
            "condition DE and R reference validation."
        ),
        parameters={"celltype_key": celltype_key, "condition_key": condition_key, "min_cells": min_cells,
                    "obs_keys": keys.record(), "design_columns": design_keys},
        evidence={"n_celltypes": len(summary_rows), "n_input_cells": int(adata.n_obs), "n_eligible_cells": len(obs),
                  "cell_filter": "annotation_status in {accepted, labeled}" if status_filtered else "none (no annotation_status)"},
    )


if __name__ == "__main__":
    main()
