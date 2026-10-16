"""Explicit assay routing and spatial feature/quality boundaries."""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np
import pandas as pd
from scipy import sparse

from .config import deep_get

SPATIAL_ASSAYS = {"xenium", "cosmx", "visium", "legacy_st"}
TARGETED_ASSAYS = {"xenium", "cosmx"}


def assay_profile(cfg: Mapping[str, Any]) -> str:
    assay = str(deep_get(cfg, "run.assay", "scrna")).lower()
    if assay not in {"scrna", *SPATIAL_ASSAYS}:
        raise ValueError(f"Unsupported run.assay: {assay}; choose scrna, xenium, cosmx, visium or legacy_st.")
    return assay


def gene_expression_mask(adata, cfg: Mapping[str, Any]) -> np.ndarray:
    """Use declared feature metadata; gene-like names do not prove gene targets."""
    key = deep_get(cfg, "qc.gene_expression_key", "is_gene_expression")
    if key in adata.var:
        values = adata.var[key]
        if values.isna().any() or not values.isin([True, False, 0, 1]).all():
            raise ValueError(f"var[{key!r}] must contain nonmissing boolean flags.")
        return values.to_numpy(bool)
    for key in ("feature_types", "feature_type"):
        if key in adata.var:
            if adata.var[key].isna().any():
                raise ValueError(f"Missing feature types in var[{key!r}].")
            return adata.var[key].astype(str).str.lower().eq("gene expression").to_numpy()
    if assay_profile(cfg) in TARGETED_ASSAYS and not deep_get(cfg, "qc.features_are_genes", False):
        raise ValueError(
            "Targeted panels need feature types/is_gene_expression, "
            "or an explicit qc.features_are_genes: true declaration."
        )
    return np.ones(adata.n_vars, dtype=bool)


def spatial_qc(adata, cfg: Mapping[str, Any]) -> pd.DataFrame:
    """Flag spatial observations without droplet-doublet or mitochondrial rules.

    Use permissive starting cutoffs, overridable for the tissue/panel. Retain
    all observations and counts, keeping exclusion decisions in ``final_use``.
    """
    assay = assay_profile(cfg)
    if assay not in SPATIAL_ASSAYS:
        raise ValueError("spatial_qc requires a spatial run.assay.")
    defaults = {
        "xenium": {"min_counts": 10, "min_genes": 5},
        "cosmx": {"min_counts": 20, "min_genes": 10},
        "visium": {"min_counts": 100, "min_genes": 50},
        "legacy_st": {"min_counts": 100, "min_genes": 50},
    }
    options = {**defaults[assay], **deep_get(cfg, "qc.spatial", {})}
    for required in ("min_counts", "min_genes"):
        if required not in options or not np.isfinite(options[required]) or options[required] < 0:
            raise ValueError(f"qc.spatial.{required} must be a finite nonnegative value.")
    sample_key = deep_get(cfg, "keys.sample", "sample_id")
    if sample_key not in adata.obs or adata.obs[sample_key].isna().any():
        raise ValueError(f"Spatial QC requires nonmissing obs[{sample_key!r}].")
    x = adata.layers[deep_get(cfg, "qc.counts_layer", "counts")]
    values = x.data if sparse.issparse(x) else np.asarray(x)
    if (
        not np.isfinite(values).all()
        or (values < 0).any()
        or not np.allclose(values, np.rint(values), rtol=0, atol=1e-6)
    ):
        raise ValueError("Spatial QC requires finite nonnegative integer counts.")
    genes = gene_expression_mask(adata, cfg)
    if not genes.any():
        raise ValueError("No Gene Expression features are available.")
    adata.var["is_gene_expression"] = genes
    total = np.asarray(x[:, genes].sum(axis=1)).ravel()
    detected = np.asarray((x[:, genes] > 0).sum(axis=1)).ravel()
    adata.obs["total_counts"] = total
    adata.obs["n_genes_by_counts"] = detected
    # A filtered gene-only matrix does not establish absence of background.
    control = np.asarray(x[:, ~genes].sum(axis=1)).ravel() if (~genes).any() else np.full(adata.n_obs, np.nan)
    adata.obs["non_gene_feature_counts"] = control
    adata.obs["non_gene_feature_fraction"] = control / np.maximum(total + control, 1)
    for flag, column in (("mt_gene", "pct_counts_mt"), ("ribo_gene", "pct_counts_ribo"), ("hb_gene", "pct_counts_hb")):
        mask = adata.var[flag].to_numpy(bool) & genes if flag in adata.var else np.zeros(adata.n_vars, bool)
        adata.obs[column] = (
            (100 * np.asarray(x[:, mask].sum(axis=1)).ravel() / np.maximum(total, 1)) if mask.any() else np.nan
        )
    coord_key = options.get("coordinate_key", "spatial")
    if coord_key not in adata.obsm:
        raise ValueError(f"Spatial QC requires obsm[{coord_key!r}].")
    coordinates = np.asarray(adata.obsm[coord_key])
    if coordinates.ndim != 2 or coordinates.shape[1] < 2:
        raise ValueError("Spatial coordinates must have at least two columns.")
    flags = {
        "low_counts": (total < options["min_counts"]) | (total == 0),
        "low_genes": (detected < options["min_genes"]) | (detected == 0),
        "invalid_coordinates": ~np.isfinite(coordinates).all(axis=1),
    }
    area_key = options.get("cell_area_key")
    if area_key:
        area = adata.obs[area_key].to_numpy(float)
        flags["invalid_area"] = ~np.isfinite(area) | (area <= 0)
        adata.obs["transcripts_per_area"] = np.divide(
            total, area, out=np.full(adata.n_obs, np.nan), where=~flags["invalid_area"]
        )
    for setting in ("in_tissue_key", "source_pass_key"):
        key = options.get(setting)
        if key:
            values = adata.obs[key]
            if values.isna().any() or not values.isin([True, False, 0, 1]).all():
                raise ValueError(f"obs[{key!r}] must contain nonmissing boolean flags.")
            flags[setting.removesuffix("_key") + "_fail"] = ~values.to_numpy(bool)
    fail = np.logical_or.reduce(list(flags.values()))
    for name, flag in flags.items():
        adata.obs[f"qc_{name}_flag"] = flag
    adata.obs["qc_pass"] = ~fail
    adata.obs["qc_class"] = np.where(fail, "fail", "pass")
    adata.obs["qc_failure_reason"] = [";".join(k for k, v in flags.items() if v[i]) for i in range(adata.n_obs)]
    previous = adata.obs.get("final_use", pd.Series(True, index=adata.obs_names))
    if previous.isna().any() or not previous.isin([True, False, 0, 1]).all():
        raise ValueError("Existing final_use must contain nonmissing boolean flags.")
    adata.obs["final_use"] = previous.to_numpy(bool) & ~fail
    adata.obs["doublet_score"] = np.nan
    adata.obs["doublet_call"] = "not_applicable_spatial"
    adata.uns["qc_assay"] = {
        "assay": assay,
        "mitochondrial_filter": "not_applied",
        "scrublet": "not_applicable",
        "feature_controls_available": bool((~genes).any()),
        "cell_area_key": area_key or "not_available",
        "coordinate_key": coord_key,
        "molecular_qv_filter": "upstream_not_recomputed_from_processed_matrix",
    }
    rows = []
    for sample, indices in adata.obs.groupby(sample_key, observed=True).indices.items():
        obs = adata.obs.iloc[indices]
        rows.append(
            {
                "sample_id": sample,
                "assay": assay,
                "n_obs": len(obs),
                "retained": int(obs["final_use"].sum()),
                "fail_fraction": float((obs.qc_class == "fail").mean()),
                "excluded_fraction": float(1 - obs.final_use.mean()),
                "median_counts": float(obs.total_counts.median()),
                "median_genes": float(obs.n_genes_by_counts.median()),
                "min_counts": options["min_counts"],
                "min_genes": options["min_genes"],
            }
        )
    return pd.DataFrame(rows)


def apply_final_use(adata) -> int:
    """Honor QC selection before fitting, including lazy AnnDataOOM row masks."""
    n_input = int(adata.n_obs)
    if "final_use" in adata.obs:
        flags = adata.obs["final_use"]
        if flags.isna().any() or not flags.isin([True, False, 0, 1]).all():
            raise ValueError("final_use must contain nonmissing boolean flags.")
        if not flags.all():
            adata._inplace_subset_obs(flags.to_numpy(bool))
    return n_input
