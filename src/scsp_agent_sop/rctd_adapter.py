"""RCTD-py 0.3.9 API adapter with explicit count and result semantics."""

from __future__ import annotations

import os
from importlib.metadata import version
from pathlib import Path
from time import perf_counter
from typing import Any, Mapping

import numpy as np
import pandas as pd
from scipy import sparse

from .assays import gene_expression_mask
from .config import deep_get
from .storage import write_json, write_table


def selected_weights(result, mode: str) -> np.ndarray:
    """Multi/doublet ``weights`` is the initial full fit, not the sparse result."""
    if mode == "full":
        return np.asarray(result.weights)
    out = np.zeros_like(result.weights)
    if mode == "multi":
        indices, weights = result.cell_type_indices, result.sub_weights
    elif mode == "doublet":
        indices = np.column_stack([result.first_type, result.second_type])
        weights = result.weights_doublet
    else:
        raise ValueError(f"Unsupported RCTD mode: {mode}")
    for row in range(len(out)):
        valid = indices[row] >= 0
        # add.at also preserves upstream results if a pair has duplicate types.
        np.add.at(out[row], indices[row, valid], weights[row, valid])
    return out


def _count_input(adata, cfg):
    layer = deep_get(cfg, "qc.counts_layer", "counts")
    if layer not in adata.layers:
        raise ValueError(f"RCTD requires explicit integer counts in layers[{layer!r}].")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("RCTD input observation and gene IDs must be unique and aligned.")
    x = adata.layers[layer]
    values = x.data if sparse.issparse(x) else np.asarray(x)
    if (
        not np.isfinite(values).all()
        or (values < 0).any()
        or not np.allclose(values, np.rint(values), rtol=0, atol=1e-6)
    ):
        raise ValueError("RCTD requires finite nonnegative integer counts, not log/transformed expression.")
    keep = adata.obs.get("final_use", pd.Series(True, index=adata.obs_names))
    if keep.isna().any() or not keep.isin([True, False, 0, 1]).all():
        raise ValueError("RCTD final_use must contain nonmissing boolean flags.")
    result = adata[keep.to_numpy(bool), gene_expression_mask(adata, cfg)].copy()
    result.X = result.layers[layer].copy()
    return result


def run_rctd_api(
    spatial_path: Path, reference_path: Path, output_path: Path, cfg: Mapping[str, Any], mode: str
) -> dict[str, Any]:
    import anndata as ad
    import rctd
    import torch

    if version("rctd-py") != "0.3.9":
        raise RuntimeError("This adapter validates RCTD-py 0.3.9; install the managed runtime.")
    options = deep_get(cfg, "spatial.rctd", {})
    started = perf_counter()
    spatial = _count_input(ad.read_h5ad(spatial_path), cfg)
    # Reference is scRNA: its own full native count totals precede intersecting
    # genes with spatial data. Never select reference HVGs to reduce input size.
    reference_cfg = {"run": {"assay": "scrna"}, "qc": {"counts_layer": deep_get(cfg, "qc.counts_layer", "counts")}}
    reference_adata = _count_input(ad.read_h5ad(reference_path), reference_cfg)
    cell_type_key = options.get("cell_type_key", "cell_type_lvl3")
    labels = reference_adata.obs[cell_type_key]
    if labels.isna().any() or labels.astype(str).str.strip().eq("").any():
        raise ValueError("RCTD reference cell types must be nonmissing reviewed identities.")
    if "annotation_status" in reference_adata.obs:
        accepted = reference_adata.obs.annotation_status.astype(str).eq("accepted")
        reference_adata = reference_adata[accepted].copy()
    reference = rctd.Reference(
        reference_adata,
        cell_type_col=cell_type_key,
        cell_min=int(options.get("reference_cell_min", 25)),
        n_max_cells=int(options.get("reference_max_cells", 3000)),
        min_UMI=int(options.get("reference_min_umi", 100)),
    )
    held_out = set(options.get("exclude_genes", []))
    if held_out:
        keep = np.array([gene not in held_out for gene in reference.gene_names])
        reference.profiles = reference.profiles[keep]
        reference.gene_names = [gene for gene, use in zip(reference.gene_names, keep) if use]
        reference.n_genes = int(keep.sum())
    q_path = options.get("q_matrices")
    previous_q = os.environ.get("RCTD_Q_MATRICES")
    if q_path:
        os.environ["RCTD_Q_MATRICES"] = str(Path(q_path).expanduser().resolve())
    device = options.get("device", "auto")
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    # Full-mode solver uses its own fixed iteration count upstream. Do not
    # advertise config fields that the selected solver does not actually use.
    config = rctd.RCTDConfig(
        device=device,
        compile=False,
        dtype="float64",
        UMI_min=int(options.get("umi_min", options.get("umi_min_low_resolution", 100))),
        MAX_MULTI_TYPES=int(options.get("max_multi_types", 4)),
    )
    try:
        model = rctd.RCTD(spatial, reference, config)
        model.fit_platform_effects(sigma_override=options.get("sigma_override"))
        arguments = dict(
            spatial_counts=model.counts,
            spatial_numi=model.nUMI,
            norm_profiles=model.norm_profiles,
            cell_type_names=reference.cell_type_names,
            q_mat=model.q_mat,
            sq_mat=model.sq_mat,
            x_vals=model.x_vals,
            batch_size=int(options.get("batch_size", 512)),
            device=device,
            pixel_mask=model._pixel_mask,
        )
        if mode != "full":
            arguments["config"] = config
        result = getattr(rctd, f"run_{mode}_mode")(**arguments)
    finally:
        if q_path:
            if previous_q is None:
                os.environ.pop("RCTD_Q_MATRICES", None)
            else:
                os.environ["RCTD_Q_MATRICES"] = previous_q
    fitted = spatial.obs_names[np.asarray(result.pixel_mask, bool)]
    raw = selected_weights(result, mode)
    if raw.shape != (len(fitted), len(result.cell_type_names)) or not np.isfinite(raw).all():
        raise ValueError("RCTD output weights are nonfinite or fail observation/type alignment.")
    table = pd.DataFrame(raw, index=fitted, columns=result.cell_type_names)
    table.index.name = "spatial_unit_id"
    positive = table.clip(lower=0)
    normalized = positive.div(positive.sum(axis=1).replace(0, np.nan), axis=0)
    diagnostics = pd.DataFrame(index=fitted)
    diagnostics["raw_weight_sum"] = table.sum(axis=1)
    diagnostics["negative_weight_mass"] = -table.clip(upper=0).sum(axis=1)
    if hasattr(result, "converged"):
        diagnostics["converged"] = result.converged
    if hasattr(result, "spot_class"):
        diagnostics["spot_class"] = [rctd.SPOT_CLASS_NAMES[int(k)] for k in result.spot_class]
    if hasattr(result, "n_types"):
        diagnostics["n_types"] = result.n_types
    folder = output_path.parent
    folder.mkdir(parents=True, exist_ok=True)
    write_table(table.reset_index(), folder / f"{mode}_weights_raw.parquet")
    write_table(normalized.reset_index(), folder / f"{mode}_weights_visualization.parquet")
    write_table(diagnostics.reset_index(), folder / f"{mode}_diagnostics.parquet")
    if mode != "full":
        write_table(
            pd.DataFrame(result.weights, index=fitted, columns=result.cell_type_names).reset_index(),
            folder / f"{mode}_weights_initial_full.parquet",
        )
    output = spatial[fitted].copy()
    output.obsm["rctd_weights_raw"] = raw
    output.obsm["rctd_proportions_visualization"] = normalized.to_numpy()
    output.uns["rctd_cell_types"] = list(result.cell_type_names)
    info = {
        "method": "rctd-py API",
        "version": version("rctd-py"),
        "mode": mode,
        "device": device,
        "compile": False,
        "dtype": "float64",
        "sigma": int(model.sigma),
        "genes_bulk": len(model._gene_list_bulk),
        "genes_reg": len(model._gene_list_reg),
        "observations_after_qc": spatial.n_obs,
        "observations_fitted": len(fitted),
        "excluded_genes": sorted(held_out),
        "reference_types": reference.cell_type_names,
        "reference_annotation_status_checked": "annotation_status" in reference_adata.obs,
        "q_matrices": q_path or previous_q or "upstream_cache_or_download",
        "elapsed_seconds": perf_counter() - started,
        "weight_semantics": "raw fitted weights; nonnegative row-normalized visualization is not absolute cell counts",
        "annotation_status": "deconvolution_estimates_not_committed_cell_identities",
    }
    output.uns["rctd_method"] = info
    output.write_h5ad(output_path)
    write_json(info, folder / f"{mode}_completion.json")
    return info
