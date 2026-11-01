"""RCTD native-count fitting, reusable profiles and fit/export recovery."""
from __future__ import annotations

import hashlib
import inspect
import json
import os
from importlib.metadata import version
from pathlib import Path
from time import perf_counter
from typing import Any, Mapping

import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits, threadpool_info

from scoop.io import atomic_json, sha256_file
from .assays import gene_expression_mask
from .config import deep_get
from .storage import write_h5ad_atomic, write_table

PRESETS = {"fast": {"N_fit": 50, "N_epoch": 4}, "standard": {"N_fit": 1000, "N_epoch": 8}}


def selected_weights(result, mode: str) -> np.ndarray:
    """Select the fitted mode's weights, including its sparse type assignment."""
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
        np.add.at(out[row], indices[row, valid], weights[row, valid])
    return out


def _attach_registry(adata, path: str | None) -> None:
    if not path:
        return
    registry = pd.read_parquet(path)
    if "source_row" in registry:
        registry = registry.sort_values("source_row")
        if not np.array_equal(registry.source_row, np.arange(adata.n_obs)):
            raise ValueError("Observation registry must cover the exact native row axis")
    if "source_axis_id" not in registry or not np.array_equal(registry.source_axis_id.astype(str), adata.obs_names.astype(str)):
        raise ValueError("Observation registry IDs/order differ from native counts")
    registry = registry.set_index(pd.Index(adata.obs_names), drop=False)
    for key in registry:
        if key in adata.obs:
            old, new = adata.obs[key].astype(object), registry[key].astype(object)
            if not ((old.eq(new)) | (old.isna() & new.isna())).all():
                raise ValueError(f"Observation registry conflicts with native metadata: {key}")
        else:
            values = registry[key]
            if (values.dtype == object or isinstance(values.dtype, pd.StringDtype)) and values.dropna().map(lambda x: isinstance(x, str)).all():
                values = pd.Categorical(values)
            adata.obs[key] = values


def _count_input(adata, cfg, options=None, role="spatial"):
    options = options or {}
    _attach_registry(adata, options.get(f"{role}_registry"))
    layer = options.get(f"{role}_counts_layer", deep_get(cfg, "qc.counts_layer", "counts"))
    if layer != "X" and layer not in adata.layers:
        raise ValueError(f"RCTD requires explicit integer counts in layers[{layer!r}], or an explicit {role}_counts_layer: X")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("RCTD observation/gene IDs must be unique; reconcile native namespaces/duplicate genes explicitly")
    x = adata.X if layer == "X" else adata.layers[layer]
    values = x.data if sparse.issparse(x) else np.asarray(x)
    if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.rint(values)).all():
        raise ValueError("RCTD requires finite nonnegative integer counts, not normalized/log expression")
    genes = gene_expression_mask(adata, cfg)
    if not genes.any():
        raise ValueError("No measured RNA features are available")
    full = adata[:, genes].copy()
    full.X = x[:, genes].copy()
    total = np.asarray(full.X.sum(axis=1)).ravel().astype(np.float64)
    exposure_key = options.get(f"{role}_exposure_key")
    exposure = full.obs[exposure_key].to_numpy(dtype=float) if exposure_key else total
    if (not np.isfinite(exposure).all() or (exposure < total).any()
            or not np.equal(exposure, np.rint(exposure)).all()):
        raise ValueError("Declared full native exposure must be finite integer depth >= measured count mass")
    if not exposure_key and "native_full_gene_UMI" in full.obs:
        declared = pd.to_numeric(full.obs.native_full_gene_UMI, errors="raise").to_numpy(dtype=float)
        measured = np.isfinite(declared)
        if not np.array_equal(declared[measured], total[measured]):
            raise ValueError("Native count depth differs from the bound observation registry")
    status = pd.DataFrame(index=full.obs_names)
    status.index.name = "spatial_unit_id"
    status["source_selected"] = True
    status["source_selection_assessed"] = False
    status["source_selection_partial"] = False
    for key in options.get("observation_mask_keys", ["final_use", "analysis_use", "source_author_retained", "coordinate_qc_finite"]):
        if key not in full.obs:
            continue
        flags = full.obs[key]
        if not flags.dropna().isin([True, False, 0, 1]).all():
            raise ValueError(f"obs[{key!r}] must be boolean; missing is unassessed, not pass")
        status["source_selected"] &= ~flags.eq(False).fillna(False).to_numpy(bool)
        status["source_selection_assessed"] |= flags.notna().to_numpy()
        status["source_selection_partial"] |= flags.isna().to_numpy()
    if "in_tissue" in full.obs:
        if not full.obs.in_tissue.dropna().isin([0, 1, True, False]).all():
            raise ValueError("in_tissue must use boolean/0/1 values")
        status["source_selected"] &= ~full.obs.in_tissue.eq(0).fillna(False).to_numpy(bool)
        status["source_selection_assessed"] |= full.obs.in_tissue.notna().to_numpy()
        status["source_selection_partial"] |= full.obs.in_tissue.isna().to_numpy()
    status["native_full_exposure"] = exposure
    return full, status


def _reference(reference_path, cfg, options):
    import anndata as ad
    import rctd

    if options.get("reference_profiles"):
        if not options.get("reference_manifest"):
            raise ValueError("reference_profiles requires its reviewed reference_manifest")
        profile_path = Path(options["reference_profiles"])
        metadata = json.loads(Path(options["reference_manifest"]).read_text())
        if metadata.get("schema_version") != "scoop.reference-profiles.v1":
            raise ValueError("Cached reference requires scoop.reference-profiles.v1 metadata")
        if (metadata.get("unit") != "mean_counts_per_full_native_exposure"
                or metadata.get("labels_reviewed") is not True
                or metadata.get("native_full_exposure_preserved") is not True):
            raise ValueError("Reference cache needs reviewed component labels and full-native count/exposure provenance")
        if sha256_file(profile_path) != metadata["profiles_sha256"]:
            raise ValueError("Cached reference profiles differ from the declared hash")
        with np.load(profile_path, allow_pickle=False) as z:
            reference = rctd.Reference.__new__(rctd.Reference)
            reference.profiles = z["profiles"].copy()
            reference.gene_names = z["gene_names"].astype(str).tolist()
            reference.cell_type_names = z["cell_type_names"].astype(str).tolist()
        reference.n_genes, reference.n_types = reference.profiles.shape
        training = metadata.get("training_cells", {})
        review = "reviewed_component_labels"
    else:
        reference_cfg = {"run": {"assay": "scrna"}, "qc": {"counts_layer": deep_get(cfg, "qc.counts_layer", "counts")}}
        full, status = _count_input(ad.read_h5ad(reference_path), reference_cfg, options, "reference")
        if not np.array_equal(status.native_full_exposure, np.asarray(full.X.sum(axis=1)).ravel()):
            raise ValueError("Projected reference counts require native-exposure profiles; use reference_profiles with its provenance manifest")
        eligible = status.source_selected.to_numpy(bool)
        mask_key = options.get("reference_mask_key")
        if mask_key:
            eligible &= full.obs[mask_key].eq(True).fillna(False).to_numpy(bool)
            review = "explicit_reference_component_mask"
        elif "annotation_status" in full.obs:
            eligible &= full.obs.annotation_status.isin(["accepted", "labeled"]).to_numpy()
            if "annotation_reviewed" in full.obs:
                eligible &= full.obs.annotation_reviewed.eq(True).fillna(False).to_numpy(bool)
            review = "reviewed_annotation_status"
        else:
            raise ValueError("Reference requires reviewed annotation_status or an explicit reviewed reference_mask_key")
        full = full[eligible].copy()
        key = options.get("cell_type_key", "cell_type_lvl3")
        labels = full.obs[key]
        if labels.isna().any() or labels.astype(str).str.strip().eq("").any():
            raise ValueError("Reference component names must be nonmissing")
        expected = set(labels.astype(str))
        reference = rctd.Reference(full, cell_type_col=key,
            cell_min=int(options.get("reference_cell_min", 1)),
            n_max_cells=int(options.get("reference_max_cells", 3000)),
            min_UMI=int(options.get("reference_min_umi", 100)))
        missing = sorted(expected - set(reference.cell_type_names))
        if missing:
            raise ValueError(f"RCTD dropped declared reference components: {missing}; qualify their counts/thresholds explicitly")
        depth = np.asarray(full.X.sum(axis=1)).ravel()
        qualified = labels.astype(str)[depth >= int(options.get("reference_min_umi", 100))]
        training = qualified.value_counts().clip(upper=int(options.get("reference_max_cells", 3000))).to_dict()
        metadata = {"coverage_complete": False, "gaps": ["whole-tissue reference coverage not assessed"]}
    if (reference.n_types < 2 or reference.n_genes < 2
            or len(set(reference.gene_names)) != reference.n_genes
            or len(reference.gene_names) != reference.n_genes
            or len(set(reference.cell_type_names)) != reference.n_types
            or len(reference.cell_type_names) != reference.n_types
            or not np.isfinite(reference.profiles).all() or (reference.profiles < 0).any()
            or (reference.profiles.sum(axis=0) <= 0).any()):
        raise ValueError("Reference profile axes, values or component exposure are invalid")
    held_out = set(options.get("exclude_genes", []))
    keep = np.array([g not in held_out for g in reference.gene_names])
    reference.profiles = reference.profiles[keep]
    reference.gene_names = [g for g, use in zip(reference.gene_names, keep) if use]
    reference.n_genes = int(keep.sum())
    if not reference.n_genes or (reference.profiles.sum(axis=0) <= 0).any():
        raise ValueError("Gene exclusion removed all expression of a declared reference component")
    return reference, {"review": review, "training_cells": training,
        "rare_components": [t for t, n in training.items() if n < 25],
        "coverage_complete": bool(metadata.get("coverage_complete", False)), "gaps": metadata.get("gaps", [])}


def _fit_payload(spatial, status, reference, options, mode):
    import rctd
    import torch

    threads = int(options.get("threads", 4))
    torch.set_num_threads(threads)
    if torch.get_num_interop_threads() != 1:
        torch.set_num_interop_threads(1)
    seed = int(options.get("seed", 271030))
    np.random.seed(seed)
    torch.manual_seed(seed)
    preset = options.get("preset", "fast")
    if preset not in PRESETS:
        raise ValueError("RCTD preset must be fast or standard")
    budget = {**PRESETS[preset], **{k: int(options[k]) for k in ["N_fit", "N_epoch"] if k in options}}
    minimum = int(options.get("umi_min", options.get("umi_min_low_resolution", 100)))
    rows = np.flatnonzero(status.source_selected.to_numpy(bool) & (status.native_full_exposure >= minimum).to_numpy())
    if not len(rows):
        raise ValueError("No spatial observations pass source selection/native exposure")
    device = options.get("device", "auto")
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    # Native exposure is selected before construction and restored before any
    # calibration, including when X is a measured common-gene count cache.
    selection = options.get("gene_selection", "de")
    if selection not in {"de", "all_shared"}:
        raise ValueError("gene_selection must be de or the explicit all_shared sensitivity")
    gene_options = {} if selection == "de" else dict(fc_cutoff=-1e9, fc_cutoff_reg=-1e9,
                                                    gene_cutoff=-1e-12, gene_cutoff_reg=-1e-12, MIN_OBS=0)
    config = rctd.RCTDConfig(device=device, compile=False, dtype="float64", UMI_min=0,
        MAX_MULTI_TYPES=int(options.get("max_multi_types", 4)), **gene_options, **budget)
    with threadpool_limits(limits=threads):
        model = rctd.RCTD(spatial[rows].copy(), reference, config)
        if not model._pixel_mask.all():
            raise ValueError("Backend unexpectedly filtered a preselected count view")
        model.nUMI = status.native_full_exposure.iloc[rows].to_numpy(dtype=float)
        model.fit_platform_effects(sigma_override=options.get("sigma_override"))
        if selection == "all_shared" and model._gene_list_reg != sorted(g for g in model.common_genes if "mt-" not in g):
            raise ValueError("The fixed shared-gene sensitivity did not retain its declared gene axis")
        args = dict(spatial_counts=model.counts, spatial_numi=model.nUMI, norm_profiles=model.norm_profiles,
            cell_type_names=reference.cell_type_names, q_mat=model.q_mat, sq_mat=model.sq_mat,
            x_vals=model.x_vals, batch_size=int(options.get("batch_size", 512)), device=device,
            pixel_mask=model._pixel_mask)
        if mode != "full":
            args["config"] = config
        result = getattr(rctd, f"run_{mode}_mode")(**args)
        actual_pools = threadpool_info()
    mask = np.asarray(result.pixel_mask, bool) if result.pixel_mask is not None else np.ones(len(rows), bool)
    fitted = rows[mask]
    raw = selected_weights(result, mode)
    if raw.shape != (len(fitted), reference.n_types) or not np.isfinite(raw).all():
        raise ValueError("RCTD weights fail finite observation/type alignment")
    if list(result.cell_type_names) != list(reference.cell_type_names):
        raise ValueError("Backend changed the reference component axis")
    arrays = {"raw": raw, "fitted_rows": fitted, "profiles": np.asarray(model.norm_profiles),
        "fit_genes": np.asarray(model._gene_list_reg, dtype=str), "types": np.asarray(reference.cell_type_names, dtype=str)}
    for attr in ["converged", "spot_class", "n_types"]:
        if hasattr(result, attr):
            arrays[attr] = np.asarray(getattr(result, attr))
    if mode != "full":
        arrays["initial_full_weights"] = np.asarray(result.weights)
    info = {"method": "rctd-py API", "version": version("rctd-py"), "mode": mode, "device": device,
        "compile": False, "dtype": "float64", "seed": seed, "threads": threads,
        "actual_torch_intra_threads": torch.get_num_threads(),
        "actual_torch_interop_threads": torch.get_num_interop_threads(), "actual_threadpools_during_fit": actual_pools,
        "preset": preset, **budget, "gene_selection": selection, "effective_native_umi_min": minimum,
        "effective_full_solver_max_iter": 50 if mode == "full" else "upstream_mode_policy",
        "ignored_requested_max_iter": options.get("max_iter", 50) if mode == "full" else "not_applicable",
        "sigma": int(model.sigma), "genes_bulk": len(model._gene_list_bulk), "genes_reg": len(model._gene_list_reg),
        "observations_input": spatial.n_obs, "observations_after_source_selection": int(status.source_selected.sum()),
        "observations_fitted": len(fitted), "reference_types": list(reference.cell_type_names),
        "native_exposure_before_gene_intersection": True, "source_masks_before_calibration": True,
        "excluded_genes": list(options.get("exclude_genes", []))}
    return arrays, info


def _lookup_path(options, rctd):
    supplied = options.get("q_matrices") or os.environ.get("RCTD_Q_MATRICES")
    if supplied:
        path = Path(supplied).expanduser().resolve()
        return path / "q_matrices.npz" if path.is_dir() else path
    path = Path(rctd.__file__).parent / "data/q_matrices.npz"
    if path.is_file():
        return path
    path = Path.home() / ".cache/rctd/q_matrices.npz"
    if not path.is_file():
        from rctd._likelihood import load_cached_q_matrices
        load_cached_q_matrices()
    return path


def run_rctd_api(spatial_path: Path, reference_path: Path | None, output_path: Path,
                 cfg: Mapping[str, Any], mode: str) -> dict[str, Any]:
    import anndata as ad
    import rctd
    import tempfile

    if version("rctd-py") != "0.4.0":
        raise RuntimeError("This adapter requires the validated rctd-py0.4.0 locked runtime")
    if mode not in {"full", "multi", "doublet"}:
        raise ValueError(f"Unsupported RCTD mode: {mode}")
    started = perf_counter()
    options = dict(deep_get(cfg, "spatial.rctd", {}))
    spatial, status = _count_input(ad.read_h5ad(spatial_path), cfg, options)
    reference, reference_info = _reference(reference_path, cfg, options)
    lookup = _lookup_path(options, rctd)
    paths = {"spatial_counts": Path(spatial_path), "likelihood_lookup": lookup}
    if options.get("reference_profiles"):
        paths["reference_profiles"] = Path(options["reference_profiles"])
    else:
        paths["reference_counts"] = Path(reference_path)
    for key in ["spatial_registry", "reference_registry"]:
        if options.get(key):
            paths[key] = Path(options[key])
    input_hashes = {key: sha256_file(path) for key, path in paths.items()}
    compute_options = {k: v for k, v in options.items() if k not in {
        "reference_manifest", "reference_profiles", "spatial_registry", "reference_registry", "q_matrices", "max_iter"}}
    code = "\n".join(inspect.getsource(f) for f in [
        selected_weights, _attach_registry, _count_input, _reference, _fit_payload, gene_expression_mask])
    identity = {"schema": "scoop.rctd-fit.v1", "inputs": input_hashes, "options": compute_options,
                "qc": deep_get(cfg, "qc", {}), "assay": deep_get(cfg, "run.assay", "scrna"), "mode": mode,
                "code_sha256": hashlib.sha256(code.encode()).hexdigest(),
                "packages": {p: version(p) for p in ["rctd-py", "torch", "numpy", "scipy"]}}
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True, default=str).encode()).hexdigest()
    folder = Path(output_path).parent
    folder.mkdir(parents=True, exist_ok=True)
    cache = folder / "cache" / mode / fingerprint
    cache.mkdir(parents=True, exist_ok=True)
    checkpoint = cache / "fit.json"
    data_path = cache / "fit.npz"
    previous_q = os.environ.get("RCTD_Q_MATRICES")
    os.environ["RCTD_Q_MATRICES"] = str(lookup)
    try:
        if checkpoint.is_file():
            saved = json.loads(checkpoint.read_text())
            if saved["fit_identity"] != fingerprint or saved["arrays_sha256"] != sha256_file(data_path):
                raise ValueError("Saved RCTD computation checkpoint is corrupt or bound to different inputs")
            with np.load(data_path, allow_pickle=False) as z:
                arrays = {k: z[k].copy() for k in z.files}
            info = saved["info"]
            reused = True
        else:
            arrays, info = _fit_payload(spatial, status, reference, options, mode)
            with tempfile.NamedTemporaryFile(dir=cache, suffix=".npz", delete=False) as handle:
                temporary = Path(handle.name)
            try:
                np.savez_compressed(temporary, **arrays)
                os.replace(temporary, data_path)
            finally:
                temporary.unlink(missing_ok=True)
            atomic_json(checkpoint, {"fit_identity": fingerprint, "identity": identity,
                        "arrays_sha256": sha256_file(data_path), "info": info})
            reused = False
    finally:
        if previous_q is None:
            os.environ.pop("RCTD_Q_MATRICES", None)
        else:
            os.environ["RCTD_Q_MATRICES"] = previous_q
    rows, raw, types = arrays["fitted_rows"], arrays["raw"], arrays["types"].tolist()
    all_raw = np.full((spatial.n_obs, len(types)), np.nan)
    all_raw[rows] = raw
    positive = np.maximum(raw, 0)
    totals = positive.sum(axis=1)
    normalized = np.full_like(all_raw, np.nan)
    normalized[rows] = np.divide(positive, totals[:, None], out=np.full_like(positive, np.nan), where=totals[:, None] > 0)
    status["fitted"] = False
    status.iloc[rows, status.columns.get_loc("fitted")] = True
    if "converged" in arrays:
        status["converged"] = pd.array([pd.NA] * len(status), dtype="boolean")
        status.iloc[rows, status.columns.get_loc("converged")] = arrays["converged"].astype(bool)
    for field in ["spot_class", "n_types"]:
        if field in arrays:
            status[field] = np.nan
            status.iloc[rows, status.columns.get_loc(field)] = arrays[field]
    status["raw_weight_sum"] = np.nan
    status.iloc[rows, status.columns.get_loc("raw_weight_sum")] = raw.sum(axis=1)
    status["negative_weight_mass"] = np.nan
    status.iloc[rows, status.columns.get_loc("negative_weight_mass")] = -np.minimum(raw, 0).sum(axis=1)
    files = []
    for name, table in [("weights_raw", pd.DataFrame(all_raw, index=status.index, columns=types)),
                        ("weights_visualization", pd.DataFrame(normalized, index=status.index, columns=types)),
                        ("diagnostics", status)]:
        path = folder / f"{mode}_{name}.parquet"
        write_table(table.reset_index(), path)
        files.append(str(path.resolve()))
    profile_path = folder / f"{mode}_fitted_profiles.parquet"
    write_table(pd.DataFrame(arrays["profiles"], index=pd.Index(arrays["fit_genes"], name="gene"), columns=types).reset_index(), profile_path)
    files.append(str(profile_path.resolve()))
    if "initial_full_weights" in arrays:
        initial = np.full_like(all_raw, np.nan)
        initial[rows] = arrays["initial_full_weights"]
        path = folder / f"{mode}_weights_initial_full.parquet"
        write_table(pd.DataFrame(initial, index=status.index, columns=types).reset_index(), path)
        files.append(str(path.resolve()))
    norms = np.linalg.norm(arrays["profiles"], axis=0)
    if not np.isfinite(norms).all() or (norms <= 0).any():
        raise ValueError("Fitted reference has nonfinite or zero-expression components")
    singular = np.linalg.svd(arrays["profiles"] / norms, compute_uv=False)
    info = {**info, "reference": reference_info, "computation_reused": reused,
            "computation_checkpoint": str(checkpoint.resolve()), "fit_identity": fingerprint,
            "source_qc_assessed_observations": int(status.source_selection_assessed.sum()),
            "unfitted_observations_retained_as_NA": int((~status.fitted).sum()),
            "profile_condition_number": float(singular[0] / singular[-1]) if singular[-1] > 0 else None,
            "profile_rank_relative1e6": int((singular > singular[0] * 1e-6).sum()),
            "fine_cell_proportion_interpretation_approved": False,
            "weight_semantics": "conditional RNA allocation, not absolute cell counts",
            "elapsed_seconds": perf_counter() - started}
    output = spatial.copy()
    output.obsm["rctd_weights_raw"] = all_raw
    output.obsm["rctd_proportions_visualization"] = normalized
    output.obs["rctd_source_selected"] = status.source_selected.to_numpy()
    output.obs["rctd_fitted"] = status.fitted.to_numpy()
    output.obs["rctd_native_full_exposure"] = status.native_full_exposure.to_numpy()
    output.uns["rctd_cell_types"] = types
    output.uns["rctd_method"] = {k: info[k] for k in ["method", "version", "mode", "sigma", "fit_identity"]}
    output.uns["rctd_method"]["details_json"] = json.dumps(info)
    # Checkpoint is already committed. A serialization failure can now recover
    # by exporting the same computation, without fitting or recalibrating again.
    write_h5ad_atomic(output, output_path)
    reread = ad.read_h5ad(output_path, backed="r")
    try:
        if reread.obs_names.astype(str).tolist() != spatial.obs_names.astype(str).tolist():
            raise ValueError("Saved RCTD observation axis changed")
    finally:
        reread.file.close()
    info["outputs"] = files + [str(Path(output_path).resolve())]
    info["inputs"] = [str(p.resolve()) for p in paths.values()]
    if options.get("reference_manifest"):
        info["inputs"].append(str(Path(options["reference_manifest"]).resolve()))
    completion = folder / f"{mode}_completion.json"
    atomic_json(completion, info)
    info["outputs"].append(str(completion.resolve()))
    return info
