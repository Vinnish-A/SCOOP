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

from scoop.io import atomic_json, default_file_mode, sha256_file
from .assays import assay_profile, gene_expression_mask
from .config import deep_get
from .measurement import var_sha256
from .reference_profiles import REVIEW_STATUSES, component_review_status, native_exposure
from .storage import write_h5ad_atomic, write_table

PRESETS = {"fast": {"N_fit": 50, "N_epoch": 4}, "standard": {"N_fit": 1000, "N_epoch": 8}}
RUN_PURPOSES = ("production", "sensitivity", "pilot")
SOURCE_MASK_KEYS = ("final_use", "analysis_use", "source_author_retained", "coordinate_qc_finite")


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
    exposure, _ = native_exposure(full.obs, total, options.get(f"{role}_exposure_key"))
    status = pd.DataFrame(index=full.obs_names)
    status.index.name = "spatial_unit_id"
    status["source_selected"] = True
    status["source_selection_assessed"] = False
    status["source_selection_partial"] = False
    for key in options.get("observation_mask_keys", SOURCE_MASK_KEYS):
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
        with np.load(Path(options["reference_profiles"]), allow_pickle=False) as z:
            reference = rctd.Reference.__new__(rctd.Reference)
            reference.profiles = z["profiles"].copy()
            reference.gene_names = z["gene_names"].astype(str).tolist()
            reference.cell_type_names = z["cell_type_names"].astype(str).tolist()
        reference.n_genes, reference.n_types = reference.profiles.shape
        metadata, training, review, components = _profile_manifest(options, reference.cell_type_names)
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
            review = "annotation_status"
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
        qualified = np.asarray(full.X.sum(axis=1)).ravel() >= int(options.get("reference_min_umi", 100))
        training = labels.astype(str)[qualified].value_counts().clip(upper=int(options.get("reference_max_cells", 3000))).to_dict()
        metadata = {"coverage_complete": False, "gaps": ["whole-tissue reference coverage not assessed"]}
        statuses, sources = component_review_status(
            full.obs[qualified], pd.Categorical(labels.astype(str)[qualified], categories=reference.cell_type_names).codes,
            reference.n_types, level=options.get("annotation_level", "fine"))
        review = f"{review}; review status from {', '.join(sources) or 'no review source (evidence_pending)'}"
        components = [{"name": str(t), "review_status": s} for t, s in zip(reference.cell_type_names, statuses)]
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
    reference.heldout_profiles = pd.DataFrame(reference.profiles[~keep], index=np.asarray(reference.gene_names)[~keep],
                                              columns=reference.cell_type_names)
    reference.profiles = reference.profiles[keep]
    reference.gene_names = [g for g, use in zip(reference.gene_names, keep) if use]
    reference.n_genes = int(keep.sum())
    if not reference.n_genes or (reference.profiles.sum(axis=0) <= 0).any():
        raise ValueError("Gene exclusion removed all expression of a declared reference component")
    return reference, {"review": review, "training_cells": training,
        "rare_components": [t for t, n in training.items() if n < 25],
        "coverage_complete": bool(metadata.get("coverage_complete", False)), "gaps": metadata.get("gaps", []),
        "components": components, "profiles_sha256": metadata.get("profiles_sha256"),
        "fit_genes_sha256": var_sha256(reference.gene_names)}


def _profile_manifest(options, names):
    """Check a cached package's manifest (v1 or v2) and return its per-component review state.

    Kept outside the fit-identity hash: the manifest gates which fits may run and how they are
    labelled; it does not change the computation.
    """
    metadata = json.loads(Path(options["reference_manifest"]).read_text())
    if (metadata.get("unit") != "mean_counts_per_full_native_exposure"
            or metadata.get("native_full_exposure_preserved") is not True):
        raise ValueError("Reference cache needs full-native count/exposure provenance")
    if sha256_file(options["reference_profiles"]) != metadata["profiles_sha256"]:
        raise ValueError("Cached reference profiles differ from the declared hash")
    schema = metadata.get("schema_version")
    # A gate_override that forbids interpretation (a declared technical pilot) leaves nothing reviewed.
    forbidden = (metadata.get("gate_override") or {}).get("biological_interpretation_allowed") is False
    if schema == "scoop.reference-profiles.v1":
        # v1 has one package-level flag: reviewed means every component is reviewed.
        reviewed = metadata.get("labels_reviewed") is True and not forbidden
        components = [{"name": n, "review_status": "reviewed" if reviewed else "evidence_pending"} for n in names]
        review = ("gate_override_forbids_interpretation" if forbidden
                  else "reviewed_component_labels" if reviewed else "unreviewed_v1_package")
        return metadata, metadata.get("training_cells", {}), review, components
    if schema != "scoop.reference-profiles.v2":
        raise ValueError("Cached reference requires scoop.reference-profiles.v1 or v2 metadata")
    declared = {c["name"]: c for c in metadata["components"]}
    if sorted(declared) != sorted(names) or len(metadata["components"]) != len(names):
        raise ValueError("Reference manifest components differ from the profile component axis")
    statuses = {c["review_status"] for c in declared.values()}
    if statuses - set(REVIEW_STATUSES):
        raise ValueError(f"Unknown component review_status: {sorted(statuses - set(REVIEW_STATUSES))}")
    components = [{**declared[n], "review_status": "evidence_pending"} if forbidden else declared[n] for n in names]
    review = "gate_override_forbids_interpretation" if forbidden else "per_component_review_status"
    return metadata, {c["name"]: c["training_cells"] for c in components}, review, components


def _compute_record(cfg, options, mode):
    """Every resolved setting the fit consumes, and the device to run on.

    The fit identity hashes this record, not the raw options: labels, gates and resources stay out,
    and explicit defaults equal implicit ones.
    """
    preset = options.get("preset", "fast")
    if preset not in PRESETS:
        raise ValueError("RCTD preset must be fast or standard")
    selection = options.get("gene_selection", "de")
    if selection not in {"de", "all_shared"}:
        raise ValueError("gene_selection must be de or the explicit all_shared sensitivity")
    device = str(options.get("device", "auto"))
    if device == "auto":
        import torch
        device = "cuda" if torch.cuda.is_available() else "cpu"
    layer = deep_get(cfg, "qc.counts_layer", "counts")
    record = {**PRESETS[preset], **{k: int(options[k]) for k in ["N_fit", "N_epoch"] if k in options},
        "gene_selection": selection, "umi_min": int(options.get("umi_min", options.get("umi_min_low_resolution", 100))),
        "exclude_genes": sorted(set(map(str, options.get("exclude_genes", [])))), "seed": int(options.get("seed", 271030)),
        "sigma_override": options.get("sigma_override"), "device_class": device.split(":")[0],
        "spatial_counts_layer": options.get("spatial_counts_layer", layer),
        "spatial_exposure_key": options.get("spatial_exposure_key"),
        "observation_mask_keys": list(options.get("observation_mask_keys", SOURCE_MASK_KEYS)),
        # The inputs of gene_expression_mask, with its defaults (its source is hashed as well).
        "gene_mask": {"assay": assay_profile(cfg), "key": deep_get(cfg, "qc.gene_expression_key", "is_gene_expression"),
                      "feature_types": sorted({str(v).lower() for v in deep_get(cfg, "qc.gene_feature_types",
                                                                                 ["gene expression", "gene"])}),
                      "features_are_genes": bool(deep_get(cfg, "qc.features_are_genes", False))}}
    if mode == "multi":
        record["max_multi_types"] = int(options.get("max_multi_types", 4))
    if not options.get("reference_profiles"):
        record["reference_h5ad"] = {"counts_layer": options.get("reference_counts_layer", layer),
            "exposure_key": options.get("reference_exposure_key"),
            "cell_type_key": options.get("cell_type_key", "cell_type_lvl3"), "mask_key": options.get("reference_mask_key"),
            "cell_min": int(options.get("reference_cell_min", 1)), "max_cells": int(options.get("reference_max_cells", 3000)),
            "min_umi": int(options.get("reference_min_umi", 100))}
    return record, device


def _fit_payload(spatial, status, reference, options, mode):
    import rctd
    import torch

    threads = int(options.get("threads", 4))
    torch.set_num_threads(threads)
    if torch.get_num_interop_threads() != 1:
        torch.set_num_interop_threads(1)
    # options carries the resolved _compute_record values.
    seed = options["seed"]
    np.random.seed(seed)
    torch.manual_seed(seed)
    preset = options.get("preset", "fast")
    budget = {k: options[k] for k in ["N_fit", "N_epoch"]}
    minimum = options["umi_min"]
    rows = np.flatnonzero(status.source_selected.to_numpy(bool) & (status.native_full_exposure >= minimum).to_numpy())
    if not len(rows):
        raise ValueError("No spatial observations pass source selection/native exposure")
    device = options["device"]
    # Native exposure is selected before construction and restored before any
    # calibration, including when X is a measured common-gene count cache.
    selection = options["gene_selection"]
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
        "excluded_genes": options["exclude_genes"]}
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


def _component_summary(weights, types, components, purpose, identifiability=None):
    """One row per component: review state and interpretation, post-fit support and identifiability.

    ``weights`` are the allocated units' visualization weights; ``identifiability`` is a table already
    bound to this reference and fit gene set.
    """
    dominant = np.bincount(weights.argmax(axis=1), minlength=len(types)) / len(weights)
    table = pd.DataFrame({"component": types,
        "review_status": [c["review_status"] for c in components],
        "biological_interpretation_allowed": [purpose != "pilot" and c["review_status"] == "reviewed" for c in components],
        "parent": [c.get("parent") for c in components], "compartment": [c.get("compartment") for c in components],
        "training_cells": [c.get("training_cells") for c in components],
        "holdout_status": [c.get("holdout_status") for c in components],
        "mean_weight": weights.mean(axis=0), "max_weight": weights.max(axis=0), "dominance_fraction": dominant,
        "frac_units_weight_ge_0.1": (weights >= 0.1).mean(axis=0),
        "effectively_absent": weights.max(axis=0) < 0.01})
    if identifiability is not None:
        keep = [c for c in ["r_true_vs_est", "dominant_recall", "identifiability_verdict", "confusable_partners"]
                if c in identifiability]
        table = table.merge(identifiability[["component", *keep]], on="component", how="left")
        if "identifiability_verdict" in table:
            table["identifiability_verdict"] = table.identifiability_verdict.fillna("not_assessed")
    return table


def _heldout_validation(spatial, status, rows, weights, heldout_profiles, seed, genes_table=None, min_counts=20):
    """F09: per held-out gene, Spearman of observed y/nUMI vs sum_k w_k * profile_gk over fitted units.

    Held-out genes get no platform-effect rescaling, so only the per-gene rank is meaningful. The
    null predicts with the weights of shuffled units.
    """
    from scipy.stats import spearmanr

    genes = [g for g in heldout_profiles.index if g in spatial.var_names]
    if not genes or len(rows) < 20:
        return pd.DataFrame(), {"status": "not_assessable" if genes else "not_applicable", "n_genes_excluded": len(heldout_profiles),
                                "n_genes_in_section": len(genes)}
    y = spatial[rows, genes].X
    y = y.toarray() if sparse.issparse(y) else np.asarray(y)
    rate = y / status.native_full_exposure.to_numpy(dtype=float)[rows, None]
    profiles = heldout_profiles.loc[genes].to_numpy()
    predicted = weights @ profiles.T
    null = weights[np.random.default_rng(seed).permutation(len(weights))] @ profiles.T
    counts = y.sum(axis=0)
    table = pd.DataFrame({"gene": genes, "section_counts": counts.astype(int),
        "spearman": [spearmanr(rate[:, j], predicted[:, j]).statistic if counts[j] >= min_counts else np.nan for j in range(len(genes))],
        "spearman_shuffled_null": [spearmanr(rate[:, j], null[:, j]).statistic if counts[j] >= min_counts else np.nan
                                   for j in range(len(genes))]})
    if genes_table is not None:
        table = table.merge(genes_table[["gene", "selection", "source_component"]], on="gene", how="left")
    scored = table.dropna(subset=["spearman"])
    summary = {"status": "assessed" if len(scored) else "not_assessable", "n_genes_excluded": len(heldout_profiles),
               "n_genes_in_section": len(genes), "n_genes_scored": len(scored), "min_section_counts": min_counts,
               "median_spearman": float(scored.spearman.median()) if len(scored) else None,
               "median_shuffled_null": float(scored.spearman_shuffled_null.median()) if len(scored) else None}
    if "selection" in scored:
        summary["by_selection"] = {sel: {"n_genes": len(g), "median_spearman": float(g.spearman.median()),
                                         "median_shuffled_null": float(g.spearman_shuffled_null.median())}
                                   for sel, g in scored.groupby("selection")}
    return table, summary


def prepare_rctd(spatial_path: Path, reference_path: Path | None, output_path: Path,
                 cfg: Mapping[str, Any], mode: str) -> dict[str, Any]:
    """Every check a fit needs before it starts; ``scoop spatial --dry-run`` stops here.

    Reads and checks the native counts and the reference, applies the run-purpose gate, checks the
    shared genes and a bound identifiability table, derives the label's artifact paths and refuses a
    completion from a different fit. Nothing is fitted or written.
    """
    import anndata as ad
    import rctd

    if version("rctd-py") != "0.4.0":
        raise RuntimeError("This adapter requires the validated rctd-py0.4.0 locked runtime")
    if mode not in {"full", "multi", "doublet"}:
        raise ValueError(f"Unsupported RCTD mode: {mode}")
    options = dict(deep_get(cfg, "spatial.rctd", {}))
    purpose = options.get("run_purpose", "production")
    if purpose not in RUN_PURPOSES:
        raise ValueError(f"spatial.rctd.run_purpose must be one of {RUN_PURPOSES}")
    merge_level = options.get("merge_level")
    if merge_level not in {None, "component", "parent", "compartment"}:
        raise ValueError("spatial.rctd.merge_level must be component, parent or compartment")
    # Published names follow output_label; only the cache is keyed by fit identity.
    label = str(options.get("output_label") or mode)
    if not label or Path(label).name != label:
        raise ValueError("spatial.rctd.output_label must be a plain file-name prefix")
    folder = Path(output_path).parent
    h5ad = folder / f"rctd_{label}.h5ad"
    if Path(output_path).name != h5ad.name:
        raise ValueError(f"The RCTD H5AD of output_label {label!r} is {h5ad.name}, not {Path(output_path).name}")
    compute, device = _compute_record(cfg, options, mode)
    spatial, status = _count_input(ad.read_h5ad(spatial_path), cfg, options)
    reference, reference_info = _reference(reference_path, cfg, options)
    pending = [c["name"] for c in reference_info["components"] if c["review_status"] == "evidence_pending"]
    if pending and purpose == "production":
        raise ValueError(f"Production RCTD fits need reviewed reference components; evidence_pending: {pending}. "
                         "Declare spatial.rctd.run_purpose: pilot or sensitivity for a non-production fit")
    if merge_level in {"parent", "compartment"} and not any(c.get(merge_level) for c in reference_info["components"]):
        raise ValueError(f"merge_level {merge_level} needs a v2 reference manifest that declares component {merge_level}s")
    shared = int(np.isin(reference.gene_names, spatial.var_names.astype(str)).sum())
    if shared < reference.n_types:
        raise ValueError(f"The spatial view shares {shared} fit genes with {reference.n_types} reference components; "
                         "check that both name genes the same way")
    identifiability = None
    if options.get("identifiability_table"):
        path = Path(options["identifiability_table"])
        identifiability = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path, sep="\t")
        bound = {"profiles_sha256": reference_info["profiles_sha256"], "reference_fit_genes_sha256": reference_info["fit_genes_sha256"]}
        if any(key not in identifiability or set(identifiability[key].astype(str)) != {str(value)} for key, value in bound.items()):
            raise ValueError("identifiability_table was not computed on these reference profiles with these excluded genes; "
                             "rerun scoop spatial identifiability")
    lookup = _lookup_path(options, rctd)
    paths = {"spatial_counts": Path(spatial_path), "likelihood_lookup": lookup}
    if options.get("reference_profiles"):
        paths["reference_profiles"] = Path(options["reference_profiles"])
    else:
        paths["reference_counts"] = Path(reference_path)
    for key in ["spatial_registry", "reference_registry"]:
        if options.get(key):
            paths[key] = Path(options[key])
    code = "\n".join(inspect.getsource(f) for f in [
        selected_weights, _attach_registry, _count_input, native_exposure, _reference, _compute_record, _fit_payload,
        gene_expression_mask])
    identity = {"schema": "scoop.rctd-fit.v2", "inputs": {key: sha256_file(path) for key, path in paths.items()},
                "compute": compute, "mode": mode, "code_sha256": hashlib.sha256(code.encode()).hexdigest(),
                "packages": {p: version(p) for p in ["rctd-py", "torch", "numpy", "scipy"]}}
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True, default=str).encode()).hexdigest()
    completion = folder / f"{label}_completion.json"
    if completion.is_file() and json.loads(completion.read_text()).get("fit_identity") != fingerprint:
        raise ValueError(f"{completion} belongs to a different RCTD fit; set a distinct spatial.rctd.output_label")
    return {"options": options, "compute": compute, "device": device, "purpose": purpose, "merge_level": merge_level,
            "label": label, "folder": folder, "h5ad": h5ad, "completion": completion, "spatial": spatial, "status": status,
            "reference": reference, "reference_info": reference_info, "pending": pending, "n_shared_genes": shared,
            "identifiability": identifiability, "lookup": lookup, "paths": paths, "identity": identity,
            "fit_identity": fingerprint}


def run_rctd_api(spatial_path: Path, reference_path: Path | None, output_path: Path,
                 cfg: Mapping[str, Any], mode: str) -> dict[str, Any]:
    """Fit RCTD after ``prepare_rctd``; ``output_path`` must be ``rctd_{output_label}.h5ad``."""
    import anndata as ad
    import tempfile
    from rctd._sigma import SIGMA_ALL

    started = perf_counter()
    plan = prepare_rctd(spatial_path, reference_path, output_path, cfg, mode)
    options, compute, spatial, status, reference, reference_info, folder, label, fingerprint = (plan[k] for k in [
        "options", "compute", "spatial", "status", "reference", "reference_info", "folder", "label", "fit_identity"])
    purpose, pending, merge_level, lookup, identity = (plan[k] for k in ["purpose", "pending", "merge_level", "lookup", "identity"])
    folder.mkdir(parents=True, exist_ok=True)
    cache = folder / "cache" / mode / fingerprint
    cache.mkdir(parents=True, exist_ok=True)
    checkpoint = cache / "fit.json"
    info_path = cache / "fit_info.json"
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
            # Earlier checkpoints carried diagnostics inline.
            info = saved["info"] if "info" in saved else json.loads(info_path.read_text())["info"]
            reused = True
        elif data_path.is_file() and info_path.is_file():
            # The arrays are the commit; a crash before fit.json leaves them reusable.
            with np.load(data_path, allow_pickle=False) as z:
                arrays = {k: z[k].copy() for k in z.files}
            diagnostics = json.loads(info_path.read_text())
            if str(arrays.get("fit_identity")) != fingerprint or diagnostics["fit_identity"] != fingerprint:
                raise ValueError("Saved RCTD arrays are bound to different inputs")
            info = diagnostics["info"]
            atomic_json(checkpoint, {"fit_identity": fingerprint, "identity": identity,
                        "arrays_sha256": sha256_file(data_path)})
            reused = True
        else:
            arrays, info = _fit_payload(spatial, status, reference, {**options, **compute, "device": plan["device"]}, mode)
            # Runtime metadata (threadpool entries) can hold callables; it must not
            # be able to fail the commit of a finished fit.
            info = json.loads(json.dumps(info, default=repr))
            atomic_json(info_path, {"fit_identity": fingerprint, "info": info})
            arrays["fit_identity"] = np.asarray(fingerprint)
            with tempfile.NamedTemporaryFile(dir=cache, suffix=".npz", delete=False) as handle:
                temporary = Path(handle.name)
            try:
                np.savez_compressed(temporary, **arrays)
                temporary.chmod(default_file_mode())  # NamedTemporaryFile creates it 0600
                os.replace(temporary, data_path)
            finally:
                temporary.unlink(missing_ok=True)
            atomic_json(checkpoint, {"fit_identity": fingerprint, "identity": identity,
                        "arrays_sha256": sha256_file(data_path)})
            reused = False
    finally:
        if previous_q is None:
            os.environ.pop("RCTD_Q_MATRICES", None)
        else:
            os.environ["RCTD_Q_MATRICES"] = previous_q
    # choose_sigma starts at 100 and moves at most 8 grid steps per epoch. A sigma
    # at that reach used every epoch, so the epoch cap, not convergence, stopped it.
    sigma_step_limited = False
    if options.get("sigma_override") is None:
        epochs = compute["N_epoch"]
        start, end = (int(np.flatnonzero(SIGMA_ALL == s)[0]) for s in (100, info["sigma"]))
        sigma_step_limited = abs(end - start) >= 8 * epochs
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
        path = folder / f"{label}_{name}.parquet"
        write_table(table.reset_index(), path)
        files.append(str(path.resolve()))
    profile_path = folder / f"{label}_fitted_profiles.parquet"
    write_table(pd.DataFrame(arrays["profiles"], index=pd.Index(arrays["fit_genes"], name="gene"), columns=types).reset_index(), profile_path)
    files.append(str(profile_path.resolve()))
    if "initial_full_weights" in arrays:
        initial = np.full_like(all_raw, np.nan)
        initial[rows] = arrays["initial_full_weights"]
        path = folder / f"{label}_weights_initial_full.parquet"
        write_table(pd.DataFrame(initial, index=status.index, columns=types).reset_index(), path)
        files.append(str(path.resolve()))
    # A fitted unit without positive raw weight has no allocation: it stays NA in the tables and
    # is left out of per-component support and held-out validation.
    allocated = np.isfinite(normalized[rows]).all(axis=1)
    fitted_weights = normalized[rows][allocated]
    summary_path = folder / f"{label}_component_summary.parquet"
    summary = _component_summary(fitted_weights, types, reference_info["components"], purpose, plan["identifiability"])
    write_table(summary, summary_path)
    files.append(str(summary_path.resolve()))
    if merge_level in {"parent", "compartment"}:
        # Exact: visualization weights are linear in the positive raw weights, so a parent's
        # weight is the sum of its components' weights.
        groups = [c.get(merge_level) or c["name"] for c in reference_info["components"]]
        merged = pd.DataFrame(normalized, index=status.index, columns=types).T.groupby(groups, sort=False).sum(min_count=1).T
        path = folder / f"{label}_weights_{merge_level}.parquet"
        write_table(merged.reset_index(), path)
        files.append(str(path.resolve()))
    genes_table = None
    if options.get("reference_manifest"):
        manifest_path = Path(options["reference_manifest"])
        declared = json.loads(manifest_path.read_text()).get("heldout_validation_genes", {}).get("path")
        if declared and (manifest_path.parent / declared).is_file():
            genes_table = pd.read_csv(manifest_path.parent / declared, sep="\t")
    heldout_table, heldout_summary = _heldout_validation(spatial, status, rows[allocated], fitted_weights,
                                                         reference.heldout_profiles, compute["seed"], genes_table)
    if len(heldout_table):
        path = folder / f"{label}_heldout_gene_validation.parquet"
        write_table(heldout_table, path)
        files.append(str(path.resolve()))
    norms = np.linalg.norm(arrays["profiles"], axis=0)
    if not np.isfinite(norms).all() or (norms <= 0).any():
        raise ValueError("Fitted reference has nonfinite or zero-expression components")
    singular = np.linalg.svd(arrays["profiles"] / norms, compute_uv=False)
    info = {**info, "reference": reference_info, "computation_reused": reused,
            "computation_checkpoint": str(checkpoint.resolve()), "fit_identity": fingerprint, "output_label": label,
            "sigma_step_limited": sigma_step_limited, "run_purpose": purpose,
            "evidence_pending_components": pending,
            "biological_interpretation_allowed": bool(summary.biological_interpretation_allowed.all()),
            "interpretation_withheld_components": summary.component[~summary.biological_interpretation_allowed].tolist(),
            "merge_level": merge_level or "component", "heldout_gene_validation": heldout_summary,
            "effectively_absent_components": summary.component[summary.effectively_absent].tolist(),
            "source_qc_assessed_observations": int(status.source_selection_assessed.sum()),
            "unfitted_observations_retained_as_NA": int((~status.fitted).sum()),
            "fitted_observations_without_allocation": int((~allocated).sum()),
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
    output.uns["rctd_method"] = {k: info[k] for k in ["method", "version", "mode", "sigma", "fit_identity", "run_purpose",
                                                      "biological_interpretation_allowed"]}
    output.uns["rctd_method"]["details_json"] = json.dumps(info)
    # Checkpoint is already committed. A serialization failure can now recover
    # by exporting the same computation, without fitting or recalibrating again.
    write_h5ad_atomic(output, plan["h5ad"])
    reread = ad.read_h5ad(plan["h5ad"], backed="r")
    try:
        if reread.obs_names.astype(str).tolist() != spatial.obs_names.astype(str).tolist():
            raise ValueError("Saved RCTD observation axis changed")
    finally:
        reread.file.close()
    info["outputs"] = files + [str(plan["h5ad"].resolve())]
    info["inputs"] = [str(p.resolve()) for p in plan["paths"].values()]
    for key in ["reference_manifest", "identifiability_table"]:
        if options.get(key):
            info["inputs"].append(str(Path(options[key]).resolve()))
    atomic_json(plan["completion"], info)
    info["outputs"].append(str(plan["completion"].resolve()))
    return info
