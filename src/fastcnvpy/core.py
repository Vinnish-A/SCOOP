from __future__ import annotations

from dataclasses import asdict, dataclass
from concurrent.futures import ThreadPoolExecutor
from time import perf_counter

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.cluster.hierarchy import cut_tree, fcluster, linkage
from scipy.spatial.distance import pdist, squareform

from .config import FastCNVConfig
from .reference import chromosome_arm_order, normalize_gene_metadata, window_chrom_arm


@dataclass(frozen=True)
class FastCNVResult:
    raw_genomic_scores: pd.DataFrame
    genomic_scores: pd.DataFrame
    cell_metadata: pd.DataFrame
    arm_cnv: pd.DataFrame
    window_metadata: pd.DataFrame
    manifest: dict


@dataclass(frozen=True)
class PooledFastCNVResult:
    per_sample: dict[str, FastCNVResult]
    cell_metadata: pd.DataFrame
    arm_cnv: pd.DataFrame
    window_metadata: pd.DataFrame
    manifest: dict


def run_fastcnv(
    counts: pd.DataFrame,
    gene_metadata: pd.DataFrame,
    *,
    obs: pd.DataFrame | None = None,
    reference_var: str | None = None,
    reference_label: str | list[str] | tuple[str, ...] | None = None,
    config: FastCNVConfig | None = None,
    sample_name: str = "sample",
    genes_to_force: list[str] | tuple[str, ...] | None = None,
    chr_arms_to_force: list[str] | tuple[str, ...] | str | None = None,
    region_to_force: tuple[str, int, int] | None = None,
    compute_clusters: bool = True,
    compute_classification: bool = True,
) -> FastCNVResult:
    """Run a Python port of fastCNV's CNVCalling/CNVPerChromosomeArm path.

    `counts` must be a gene-by-cell matrix, matching Seurat's assay counts orientation.
    """

    cfg = config or FastCNVConfig()
    cfg.validate()
    timings: dict[str, float] = {}
    total_start = perf_counter()

    start = perf_counter()
    gene_names = counts.index.astype(str)
    cell_names = counts.columns.astype(str)
    if not counts.index.equals(gene_names) or not counts.columns.equals(cell_names):
        counts = counts.copy(deep=False)
        counts.index = gene_names
        counts.columns = cell_names
    n_input_genes = int(counts.shape[0])
    obs = _align_obs(obs, counts.columns)
    gmeta = normalize_gene_metadata(gene_metadata)
    common_genes = pd.Index(counts.index).intersection(pd.Index(gmeta["hgnc_symbol"]))
    if len(common_genes) == 0:
        raise ValueError("no overlapping genes between counts and fastCNV gene metadata")
    counts = counts.loc[common_genes].astype("float64", copy=False)
    reference_cells, scale_on_reference = _reference_cells(
        obs=obs,
        cell_names=pd.Index(counts.columns),
        reference_var=reference_var,
        reference_label=reference_label,
        scale_on_reference=cfg.scale_on_reference_label,
    )
    timings["prepare_seconds"] = perf_counter() - start

    start = perf_counter()
    average_expression = _average_expression(counts, reference_cells)
    selected_genes = _select_genes_like_fastcnv(
        common_genes=list(common_genes),
        average_expression=average_expression,
        gene_metadata=gmeta,
        top_n_genes=min(cfg.top_n_genes, len(common_genes)),
        genes_to_force=genes_to_force,
        chr_arms_to_force=chr_arms_to_force,
        region_to_force=region_to_force,
    )
    counts = counts.loc[selected_genes.final_selected_genes]
    timings["select_genes_seconds"] = perf_counter() - start

    start = perf_counter()
    norm = np.log2(1.0 + counts.to_numpy(dtype=np.float64, copy=False))
    norm = norm - np.nanmean(norm, axis=0, keepdims=True)
    scale_factor = _scale_factor(norm, counts.columns, reference_cells, scale_on_reference)
    norm = norm - scale_factor[:, None]
    norm = np.clip(norm, -3.0, 3.0)
    timings["normalize_seconds"] = perf_counter() - start

    start = perf_counter()
    windows = _build_genomic_windows(
        selected_genes.top_gene_metadata,
        window_size=cfg.window_size,
        window_step=cfg.window_step,
    )
    raw_scores = _compute_genomic_scores_from_array(norm, counts.index, windows)
    trimmed = _threshold_scores_like_fastcnv(
        raw_scores=raw_scores,
        cell_names=pd.Index(counts.columns),
        obs=obs,
        reference_var=reference_var,
        reference_cells=reference_cells,
        scale_on_reference=scale_on_reference,
        threshold_percentile=cfg.threshold_percentile,
    )
    timings["windows_and_threshold_seconds"] = perf_counter() - start

    start = perf_counter()
    raw_df = pd.DataFrame(raw_scores.T, index=list(windows.keys()), columns=counts.columns)
    genomic_df = pd.DataFrame(trimmed, index=list(windows.keys()), columns=counts.columns)
    cnv_fraction = (np.abs(genomic_df.to_numpy()) > 0).mean(axis=0)
    cell_meta = pd.DataFrame(index=counts.columns)
    cell_meta["cnv_fraction"] = cnv_fraction
    arm_cnv = cnv_per_chromosome_arm(genomic_df)
    cell_meta = cell_meta.join(arm_cnv.T)
    if compute_classification:
        class_df = classify_arm_cnv(arm_cnv, cfg.classification_peaks)
        cell_meta = cell_meta.join(class_df.T)
    if compute_clusters:
        clusters = cluster_cnv(genomic_df, k=cfg.cluster_k, h=cfg.cluster_h)
        cell_meta["cnv_clusters"] = clusters.reindex(cell_meta.index).astype("Int64")
        if cfg.merge_cnv:
            cell_meta["cnv_clusters"] = merge_cnv_clusters(cell_meta, merge_threshold=cfg.merge_threshold)
    window_meta = _window_metadata(windows)
    timings["summarize_seconds"] = perf_counter() - start

    timings["total_seconds"] = perf_counter() - total_start
    manifest = {
        "schema_version": "fastcnvpy.result.v1",
        "sample_name": sample_name,
        "engine": "fastcnvpy",
        "source_algorithm": "must-bioinfo/fastCNV CNVCalling-compatible port",
        "config": asdict(cfg),
        "n_input_genes": n_input_genes,
        "n_cells": int(len(counts.columns)),
        "n_common_genes": int(len(common_genes)),
        "n_final_selected_genes": int(len(selected_genes.final_selected_genes)),
        "n_top_expr_genes": int(len(selected_genes.top_expr_genes)),
        "n_windows": int(len(windows)),
        "reference_var": reference_var,
        "reference_label": list(reference_label) if isinstance(reference_label, (list, tuple)) else reference_label,
        "scale_on_reference_label": bool(scale_on_reference),
        "reference_cells": int(sum(len(v) for v in reference_cells.values())) if isinstance(reference_cells, dict) else int(len(reference_cells)),
        "timings": {key: round(value, 6) for key, value in timings.items()},
    }
    return FastCNVResult(
        raw_genomic_scores=raw_df,
        genomic_scores=genomic_df,
        cell_metadata=cell_meta,
        arm_cnv=arm_cnv,
        window_metadata=window_meta,
        manifest=manifest,
    )


def run_fastcnv_anndata(
    adata,
    gene_metadata: pd.DataFrame,
    *,
    layer: str | None = None,
    obs: pd.DataFrame | None = None,
    reference_var: str | None = None,
    reference_label: str | list[str] | tuple[str, ...] | None = None,
    config: FastCNVConfig | None = None,
    sample_name: str = "sample",
    genes_to_force: list[str] | tuple[str, ...] | None = None,
    chr_arms_to_force: list[str] | tuple[str, ...] | str | None = None,
    region_to_force: tuple[str, int, int] | None = None,
    compute_clusters: bool = True,
    compute_classification: bool = True,
    densify_all: bool = False,
) -> FastCNVResult:
    """Run fastCNV directly from an AnnData object without densifying all genes."""

    cfg = config or FastCNVConfig()
    cfg.validate()
    timings: dict[str, float] = {}
    total_start = perf_counter()

    start = perf_counter()
    matrix = adata.layers[layer] if layer else (adata.layers["counts"] if "counts" in adata.layers else adata.X)
    if densify_all and sp.issparse(matrix):
        matrix = matrix.toarray()
    cell_names = pd.Index(adata.obs_names.astype(str))
    gene_names = pd.Index(adata.var_names.astype(str))
    n_input_genes = int(len(gene_names))
    obs = _align_obs(obs if obs is not None else adata.obs.copy(), cell_names)
    gmeta = normalize_gene_metadata(gene_metadata)
    common_genes = _ordered_unique_intersection(gene_names, pd.Index(gmeta["hgnc_symbol"]))
    if len(common_genes) == 0:
        raise ValueError("no overlapping genes between counts and fastCNV gene metadata")
    first_gene_pos = _first_positions(gene_names)
    common_gene_pos = np.asarray([first_gene_pos[gene] for gene in common_genes], dtype=np.intp)
    reference_cells, scale_on_reference = _reference_cells(
        obs=obs,
        cell_names=cell_names,
        reference_var=reference_var,
        reference_label=reference_label,
        scale_on_reference=cfg.scale_on_reference_label,
    )
    timings["prepare_seconds"] = perf_counter() - start

    start = perf_counter()
    reference_pos = _cell_positions(cell_names, reference_cells)
    average_values = _average_expression_anndata_matrix(matrix, common_gene_pos, reference_pos)
    average_expression = pd.Series(average_values, index=common_genes)
    selected_genes = _select_genes_like_fastcnv(
        common_genes=list(common_genes),
        average_expression=average_expression,
        gene_metadata=gmeta,
        top_n_genes=min(cfg.top_n_genes, len(common_genes)),
        genes_to_force=genes_to_force,
        chr_arms_to_force=chr_arms_to_force,
        region_to_force=region_to_force,
    )
    final_gene_pos = np.asarray([first_gene_pos[gene] for gene in selected_genes.final_selected_genes], dtype=np.intp)
    timings["select_genes_seconds"] = perf_counter() - start

    start = perf_counter()
    counts_values = _extract_gene_by_cell_array(matrix, final_gene_pos)
    norm = np.log2(1.0 + counts_values)
    norm = norm - np.nanmean(norm, axis=0, keepdims=True)
    scale_factor = _scale_factor(norm, cell_names, reference_cells, scale_on_reference)
    norm = norm - scale_factor[:, None]
    norm = np.clip(norm, -3.0, 3.0)
    gene_index = pd.Index(selected_genes.final_selected_genes)
    timings["normalize_seconds"] = perf_counter() - start

    start = perf_counter()
    windows = _build_genomic_windows(
        selected_genes.top_gene_metadata,
        window_size=cfg.window_size,
        window_step=cfg.window_step,
    )
    raw_scores = _compute_genomic_scores_from_array(norm, gene_index, windows)
    trimmed = _threshold_scores_like_fastcnv(
        raw_scores=raw_scores,
        cell_names=cell_names,
        obs=obs,
        reference_var=reference_var,
        reference_cells=reference_cells,
        scale_on_reference=scale_on_reference,
        threshold_percentile=cfg.threshold_percentile,
    )
    timings["windows_and_threshold_seconds"] = perf_counter() - start

    result = _summarize_fastcnv_result(
        raw_scores=raw_scores,
        trimmed=trimmed,
        windows=windows,
        cell_names=cell_names,
        compute_clusters=compute_clusters,
        compute_classification=compute_classification,
        cfg=cfg,
    )
    timings["summarize_seconds"] = result.manifest.pop("_summarize_seconds")
    timings["total_seconds"] = perf_counter() - total_start
    result.manifest.update(
        {
            "schema_version": "fastcnvpy.result.v1",
            "sample_name": sample_name,
            "engine": "fastcnvpy",
            "source_algorithm": "must-bioinfo/fastCNV CNVCalling-compatible AnnData port",
            "config": asdict(cfg),
            "n_input_genes": n_input_genes,
            "n_cells": int(len(cell_names)),
            "n_common_genes": int(len(common_genes)),
            "n_final_selected_genes": int(len(selected_genes.final_selected_genes)),
            "n_top_expr_genes": int(len(selected_genes.top_expr_genes)),
            "n_windows": int(len(windows)),
            "reference_var": reference_var,
            "reference_label": list(reference_label) if isinstance(reference_label, (list, tuple)) else reference_label,
            "scale_on_reference_label": bool(scale_on_reference),
            "reference_cells": int(sum(len(v) for v in reference_cells.values())) if isinstance(reference_cells, dict) else int(len(reference_cells)),
            "timings": {key: round(value, 6) for key, value in timings.items()},
        }
    )
    return result


def run_fastcnv_pooled_anndata(
    adata,
    gene_metadata: pd.DataFrame,
    *,
    sample_key: str = "sample_id",
    layer: str | None = None,
    obs: pd.DataFrame | None = None,
    reference_var: str | None = None,
    reference_label: str | list[str] | tuple[str, ...] | None = None,
    config: FastCNVConfig | None = None,
    sample_name: str = "pooled",
    genes_to_force: list[str] | tuple[str, ...] | None = None,
    chr_arms_to_force: list[str] | tuple[str, ...] | str | None = None,
    region_to_force: tuple[str, int, int] | None = None,
    compute_clusters: bool = True,
    compute_classification: bool = True,
    densify_all: bool = False,
    n_jobs: int = 1,
    min_reference_cells_per_sample: int = 5,
) -> PooledFastCNVResult:
    """Run fastCNV's list-style pooled-reference workflow on a merged AnnData.

    This mirrors fastCNV's `CNVCallingList`: samples share common genes, selected
    genes, scale factor, and quantile thresholds, while genomic scores are
    computed and stored per sample.
    """

    cfg = config or FastCNVConfig()
    cfg.validate()
    timings: dict[str, float] = {}
    total_start = perf_counter()

    start = perf_counter()
    matrix = adata.layers[layer] if layer else (adata.layers["counts"] if "counts" in adata.layers else adata.X)
    if densify_all and sp.issparse(matrix):
        matrix = matrix.toarray()
    cell_names = pd.Index(adata.obs_names.astype(str))
    gene_names = pd.Index(adata.var_names.astype(str))
    n_input_genes = int(len(gene_names))
    obs = _align_obs(obs if obs is not None else adata.obs.copy(), cell_names)
    if sample_key not in obs.columns:
        raise ValueError(f"sample_key {sample_key!r} is not present in obs")
    gmeta = normalize_gene_metadata(gene_metadata)
    common_genes = _ordered_unique_intersection(gene_names, pd.Index(gmeta["hgnc_symbol"]))
    if len(common_genes) == 0:
        raise ValueError("no overlapping genes between counts and fastCNV gene metadata")
    first_gene_pos = _first_positions(gene_names)
    common_gene_pos = np.asarray([first_gene_pos[gene] for gene in common_genes], dtype=np.intp)
    sample_groups = {
        str(sample): cell_names.get_indexer(pd.Index(names).astype(str)).astype(np.intp)
        for sample, names in obs.groupby(sample_key, sort=True).groups.items()
    }
    sample_groups = {sample: idx[idx >= 0] for sample, idx in sample_groups.items() if len(idx[idx >= 0]) > 0}
    reference_by_sample, reference_by_label_sample, scale_on_reference = _reference_cells_by_sample(
        obs=obs,
        cell_names=cell_names,
        sample_groups=sample_groups,
        reference_var=reference_var,
        reference_label=reference_label,
        scale_on_reference=cfg.scale_on_reference_label,
        min_cells=min_reference_cells_per_sample,
    )
    timings["prepare_seconds"] = perf_counter() - start

    start = perf_counter()
    average_expression = _pooled_average_expression(
        matrix=matrix,
        common_gene_pos=common_gene_pos,
        common_genes=common_genes,
        sample_groups=sample_groups,
        reference_by_sample=reference_by_sample,
        scale_on_reference=scale_on_reference,
    )
    selected_genes = _select_genes_like_fastcnv(
        common_genes=list(common_genes),
        average_expression=average_expression,
        gene_metadata=gmeta,
        top_n_genes=min(cfg.top_n_genes, len(common_genes)),
        genes_to_force=genes_to_force,
        chr_arms_to_force=chr_arms_to_force,
        region_to_force=region_to_force,
    )
    final_gene_pos = np.asarray([first_gene_pos[gene] for gene in selected_genes.final_selected_genes], dtype=np.intp)
    gene_index = pd.Index(selected_genes.final_selected_genes)
    timings["select_genes_seconds"] = perf_counter() - start

    start = perf_counter()
    scale_factor = _pooled_scale_factor(
        matrix=matrix,
        final_gene_pos=final_gene_pos,
        sample_groups=sample_groups,
        reference_by_sample=reference_by_sample,
        reference_by_label_sample=reference_by_label_sample,
        scale_on_reference=scale_on_reference,
    )
    timings["scale_factor_seconds"] = perf_counter() - start

    start = perf_counter()
    windows = _build_genomic_windows(
        selected_genes.top_gene_metadata,
        window_size=cfg.window_size,
        window_step=cfg.window_step,
    )
    timings["windows_seconds"] = perf_counter() - start

    def compute_one(item: tuple[str, np.ndarray]) -> tuple[str, np.ndarray, pd.Index]:
        sample, cell_pos = item
        values = _extract_gene_by_cell_array_for_cells(matrix, cell_pos, final_gene_pos)
        norm = np.log2(1.0 + values)
        norm = norm - np.nanmean(norm, axis=0, keepdims=True)
        norm = norm - scale_factor[:, None]
        norm = np.clip(norm, -3.0, 3.0)
        raw_scores = _compute_genomic_scores_from_array(norm, gene_index, windows)
        return sample, raw_scores, cell_names.take(cell_pos)

    start = perf_counter()
    items = list(sample_groups.items())
    if n_jobs > 1 and len(items) > 1:
        with ThreadPoolExecutor(max_workers=int(n_jobs)) as pool:
            computed = list(pool.map(compute_one, items))
    else:
        computed = [compute_one(item) for item in items]
    raw_by_sample = {sample: raw for sample, raw, _ in computed}
    cells_by_sample = {sample: names for sample, _, names in computed}
    timings["sample_scores_seconds"] = perf_counter() - start

    start = perf_counter()
    threshold_source = _pooled_threshold_source(raw_by_sample, reference_by_sample, sample_groups, scale_on_reference)
    q = np.nanquantile(threshold_source, [cfg.threshold_percentile, 1.0 - cfg.threshold_percentile], axis=0)
    per_sample: dict[str, FastCNVResult] = {}
    combined_cell_meta = []
    combined_arm_cnv = []
    for sample in sorted(raw_by_sample):
        raw_scores = raw_by_sample[sample]
        trimmed_for_threshold = raw_scores.copy()
        inside = (trimmed_for_threshold >= q[0][None, :]) & (trimmed_for_threshold <= q[1][None, :])
        trimmed_for_threshold[inside] = 0.0
        result = _summarize_fastcnv_result(
            raw_scores=raw_scores,
            trimmed=trimmed_for_threshold.T,
            windows=windows,
            cell_names=cells_by_sample[sample],
            compute_clusters=compute_clusters,
            compute_classification=compute_classification,
            cfg=cfg,
        )
        result.manifest.update({"sample_name": sample, "n_cells": int(raw_scores.shape[0])})
        per_sample[sample] = result
        cell_meta = result.cell_metadata.copy()
        cell_meta[sample_key] = sample
        combined_cell_meta.append(cell_meta)
        arm = result.arm_cnv.copy()
        arm.columns = pd.MultiIndex.from_product([[sample], arm.columns], names=[sample_key, "cell"])
        combined_arm_cnv.append(arm)
    cell_metadata = pd.concat(combined_cell_meta, axis=0).reindex(cell_names)
    arm_cnv = pd.concat(combined_arm_cnv, axis=1) if combined_arm_cnv else pd.DataFrame()
    window_meta = _window_metadata(windows)
    timings["threshold_and_summarize_seconds"] = perf_counter() - start

    timings["total_seconds"] = perf_counter() - total_start
    manifest = {
        "schema_version": "fastcnvpy.pooled_result.v1",
        "sample_name": sample_name,
        "engine": "fastcnvpy",
        "source_algorithm": "must-bioinfo/fastCNV CNVCallingList-compatible pooled-reference AnnData port",
        "config": asdict(cfg),
        "sample_key": sample_key,
        "samples": sorted(sample_groups),
        "n_samples": int(len(sample_groups)),
        "n_input_genes": n_input_genes,
        "n_cells": int(len(cell_names)),
        "n_common_genes": int(len(common_genes)),
        "n_final_selected_genes": int(len(selected_genes.final_selected_genes)),
        "n_top_expr_genes": int(len(selected_genes.top_expr_genes)),
        "n_windows": int(len(windows)),
        "reference_var": reference_var,
        "reference_label": list(reference_label) if isinstance(reference_label, (list, tuple)) else reference_label,
        "scale_on_reference_label": bool(scale_on_reference),
        "reference_cells_by_sample": {sample: int(len(pos)) for sample, pos in reference_by_sample.items()},
        "n_jobs": int(n_jobs),
        "h5ad_mode": "dense" if densify_all else "sparse",
        "timings": {key: round(value, 6) for key, value in timings.items()},
    }
    return PooledFastCNVResult(
        per_sample=per_sample,
        cell_metadata=cell_metadata,
        arm_cnv=arm_cnv,
        window_metadata=window_meta,
        manifest=manifest,
    )


def _summarize_fastcnv_result(
    *,
    raw_scores: np.ndarray,
    trimmed: np.ndarray,
    windows: dict[str, list[str]],
    cell_names: pd.Index,
    compute_clusters: bool,
    compute_classification: bool,
    cfg: FastCNVConfig,
) -> FastCNVResult:
    start = perf_counter()
    raw_df = pd.DataFrame(raw_scores.T, index=list(windows.keys()), columns=cell_names)
    genomic_df = pd.DataFrame(trimmed, index=list(windows.keys()), columns=cell_names)
    cnv_fraction = (np.abs(genomic_df.to_numpy()) > 0).mean(axis=0)
    cell_meta = pd.DataFrame(index=cell_names)
    cell_meta["cnv_fraction"] = cnv_fraction
    arm_cnv = cnv_per_chromosome_arm(genomic_df)
    cell_meta = cell_meta.join(arm_cnv.T)
    if compute_classification:
        class_df = classify_arm_cnv(arm_cnv, cfg.classification_peaks)
        cell_meta = cell_meta.join(class_df.T)
    if compute_clusters:
        clusters = cluster_cnv(genomic_df, k=cfg.cluster_k, h=cfg.cluster_h)
        cell_meta["cnv_clusters"] = clusters.reindex(cell_meta.index).astype("Int64")
        if cfg.merge_cnv:
            cell_meta["cnv_clusters"] = merge_cnv_clusters(cell_meta, merge_threshold=cfg.merge_threshold)
    window_meta = _window_metadata(windows)
    return FastCNVResult(
        raw_genomic_scores=raw_df,
        genomic_scores=genomic_df,
        cell_metadata=cell_meta,
        arm_cnv=arm_cnv,
        window_metadata=window_meta,
        manifest={"_summarize_seconds": perf_counter() - start},
    )


def _ordered_unique_intersection(left: pd.Index, right: pd.Index) -> pd.Index:
    right_values = set(right.astype(str))
    seen = set()
    out = []
    for value in left.astype(str):
        if value in right_values and value not in seen:
            seen.add(value)
            out.append(value)
    return pd.Index(out)


def _first_positions(index: pd.Index) -> dict[str, int]:
    positions: dict[str, int] = {}
    for pos, value in enumerate(index.astype(str)):
        positions.setdefault(value, pos)
    return positions


def _cell_positions(cell_names: pd.Index, reference_cells) -> np.ndarray:
    if isinstance(reference_cells, dict):
        cells = [cell for values in reference_cells.values() for cell in values]
    else:
        cells = list(reference_cells)
    indexer = cell_names.get_indexer(cells)
    indexer = indexer[indexer >= 0]
    if len(indexer) == 0:
        return np.arange(len(cell_names), dtype=np.intp)
    return indexer.astype(np.intp, copy=False)


def _reference_cells_by_sample(
    *,
    obs: pd.DataFrame,
    cell_names: pd.Index,
    sample_groups: dict[str, np.ndarray],
    reference_var: str | None,
    reference_label,
    scale_on_reference: bool,
    min_cells: int,
) -> tuple[dict[str, np.ndarray], dict[str, dict[str, np.ndarray]], bool]:
    if obs is None or reference_var is None or reference_label is None or reference_var not in obs.columns:
        return {}, {}, False
    labels = [reference_label] if isinstance(reference_label, str) else list(reference_label)
    pos = {cell: i for i, cell in enumerate(cell_names)}
    out: dict[str, np.ndarray] = {}
    by_label: dict[str, dict[str, np.ndarray]] = {}
    for sample, sample_pos in sample_groups.items():
        names = cell_names.take(sample_pos)
        sub = obs.loc[names]
        sample_refs = []
        for label in labels:
            label_str = str(label)
            ref_names = pd.Index(sub.index[sub[reference_var].astype(str) == label_str])
            idx = np.asarray([pos[name] for name in ref_names if name in pos], dtype=np.intp)
            if len(idx) >= min_cells:
                by_label.setdefault(label_str, {})[sample] = idx
                sample_refs.append(idx)
        if sample_refs:
            out[sample] = np.unique(np.concatenate(sample_refs)).astype(np.intp, copy=False)
    if not out:
        return {}, {}, False
    return out, by_label, scale_on_reference


def _pooled_average_expression(
    *,
    matrix,
    common_gene_pos: np.ndarray,
    common_genes: pd.Index,
    sample_groups: dict[str, np.ndarray],
    reference_by_sample: dict[str, np.ndarray],
    scale_on_reference: bool,
) -> pd.Series:
    # fastCNV's list branch uses whole-sample means for samples that contain
    # reference cells; reference cells are used later for scaling/thresholding.
    sample_positions = {
        sample: sample_groups[sample]
        for sample in (reference_by_sample if scale_on_reference and reference_by_sample else sample_groups)
        if sample in sample_groups
    }
    averages = [
        _average_expression_anndata_matrix(matrix, common_gene_pos, pos)
        for sample, pos in sample_positions.items()
        if sample in sample_groups and len(pos) > 0
    ]
    if not averages:
        averages = [_average_expression_anndata_matrix(matrix, common_gene_pos, pos) for pos in sample_groups.values()]
    return pd.Series(np.nanmean(np.vstack(averages), axis=0), index=common_genes)


def _pooled_scale_factor(
    *,
    matrix,
    final_gene_pos: np.ndarray,
    sample_groups: dict[str, np.ndarray],
    reference_by_sample: dict[str, np.ndarray],
    reference_by_label_sample: dict[str, dict[str, np.ndarray]],
    scale_on_reference: bool,
) -> np.ndarray:
    if scale_on_reference and len(reference_by_label_sample) > 1:
        label_factors = [
            _pooled_scale_factor(
                matrix=matrix,
                final_gene_pos=final_gene_pos,
                sample_groups=sample_groups,
                reference_by_sample=label_refs,
                reference_by_label_sample={},
                scale_on_reference=True,
            )
            for label_refs in reference_by_label_sample.values()
            if label_refs
        ]
        if label_factors:
            return np.nanmedian(np.vstack(label_factors), axis=0)

    sample_factors = []
    label_positions = reference_by_sample if scale_on_reference and reference_by_sample else sample_groups
    for sample, sample_pos in sample_groups.items():
        values = _extract_gene_by_cell_array_for_cells(matrix, sample_pos, final_gene_pos)
        norm = np.log2(1.0 + values)
        norm = norm - np.nanmean(norm, axis=0, keepdims=True)
        use_pos = label_positions.get(sample, sample_pos)
        local = _positions_within_sample(sample_pos, use_pos)
        if len(local) > 0:
            sample_factors.append(np.nanmean(norm[:, local], axis=1))
    if not sample_factors:
        return np.zeros(len(final_gene_pos), dtype=np.float64)
    return np.nanmean(np.vstack(sample_factors), axis=0)


def _pooled_threshold_source(
    raw_by_sample: dict[str, np.ndarray],
    reference_by_sample: dict[str, np.ndarray],
    sample_groups: dict[str, np.ndarray],
    scale_on_reference: bool,
) -> np.ndarray:
    pieces = []
    if scale_on_reference and reference_by_sample:
        for sample, global_ref_pos in reference_by_sample.items():
            if sample not in raw_by_sample or sample not in sample_groups:
                continue
            local = _positions_within_sample(sample_groups[sample], global_ref_pos)
            if len(local) > 0:
                pieces.append(raw_by_sample[sample][local, :])
    else:
        pieces = list(raw_by_sample.values())
    if not pieces:
        pieces = list(raw_by_sample.values())
    return np.vstack(pieces)


def _positions_within_sample(sample_pos: np.ndarray, global_pos: np.ndarray) -> np.ndarray:
    local = {int(pos): i for i, pos in enumerate(sample_pos)}
    return np.asarray([local[int(pos)] for pos in global_pos if int(pos) in local], dtype=np.intp)


def _average_expression_anndata_matrix(matrix, gene_positions: np.ndarray, cell_positions: np.ndarray) -> np.ndarray:
    if sp.issparse(matrix):
        if len(cell_positions) == matrix.shape[0] and np.array_equal(cell_positions, np.arange(matrix.shape[0])):
            sub = matrix[:, gene_positions]
        else:
            sub = matrix[cell_positions, :][:, gene_positions]
        return np.asarray(sub.sum(axis=0)).ravel() / float(len(cell_positions))
    if len(cell_positions) == matrix.shape[0] and np.array_equal(cell_positions, np.arange(matrix.shape[0])):
        return np.asarray(matrix).mean(axis=0)[gene_positions]
    else:
        sub = np.asarray(matrix[cell_positions, :])
        return sub.mean(axis=0)[gene_positions]


def _extract_gene_by_cell_array(matrix, gene_positions: np.ndarray) -> np.ndarray:
    if sp.issparse(matrix):
        return matrix[:, gene_positions].T.toarray().astype(np.float64, copy=False)
    return np.asarray(matrix[:, gene_positions], dtype=np.float64).T


def _extract_gene_by_cell_array_for_cells(matrix, cell_positions: np.ndarray, gene_positions: np.ndarray) -> np.ndarray:
    if sp.issparse(matrix):
        return matrix[cell_positions, :][:, gene_positions].T.toarray().astype(np.float64, copy=False)
    return np.asarray(matrix[np.ix_(cell_positions, gene_positions)], dtype=np.float64).T


@dataclass(frozen=True)
class SelectedGenes:
    top_expr_genes: list[str]
    final_selected_genes: list[str]
    top_gene_metadata: pd.DataFrame


def _align_obs(obs: pd.DataFrame | None, cell_names: pd.Index) -> pd.DataFrame | None:
    if obs is None:
        return None
    out = obs.copy()
    out.index = out.index.astype(str)
    return out.reindex(cell_names)


def _reference_cells(
    *,
    obs: pd.DataFrame | None,
    cell_names: pd.Index,
    reference_var: str | None,
    reference_label,
    scale_on_reference: bool,
):
    if obs is None or reference_var is None or reference_label is None or reference_var not in obs.columns:
        return list(cell_names), False
    labels = [reference_label] if isinstance(reference_label, str) else list(reference_label)
    if len(labels) == 1:
        cells = list(obs.index[obs[reference_var].astype(str) == str(labels[0])])
        if len(cells) == 0:
            return list(cell_names), False
        return cells, scale_on_reference
    out = {}
    for label in labels:
        cells = list(obs.index[obs[reference_var].astype(str) == str(label)])
        if len(cells) >= 5:
            out[str(label)] = cells
    if not out:
        return list(cell_names), False
    return out, scale_on_reference


def _average_expression(counts: pd.DataFrame, reference_cells) -> pd.Series:
    if isinstance(reference_cells, dict):
        cells = [cell for values in reference_cells.values() for cell in values]
    else:
        cells = list(reference_cells)
    columns = pd.Index(counts.columns)
    idx = columns.get_indexer(cells)
    idx = idx[idx >= 0]
    values = counts.to_numpy(dtype=np.float64, copy=False)
    if len(idx) > 0:
        return pd.Series(values[:, idx].mean(axis=1), index=counts.index)
    return pd.Series(values.mean(axis=1), index=counts.index)


def _select_genes_like_fastcnv(
    *,
    common_genes: list[str],
    average_expression: pd.Series,
    gene_metadata: pd.DataFrame,
    top_n_genes: int,
    genes_to_force,
    chr_arms_to_force,
    region_to_force,
) -> SelectedGenes:
    ave = average_expression.reindex(common_genes)
    top_expr = list(ave.sort_values(ascending=False, kind="mergesort").index[:top_n_genes])
    if genes_to_force is not None:
        top_expr = _union_preserve(top_expr, [gene for gene in genes_to_force if gene in common_genes])
    if region_to_force is not None:
        chrom, start, end = region_to_force
        region = gene_metadata[
            (gene_metadata["chromosome_num"].astype(str).replace({"23": "X"}) == str(chrom))
            & (gene_metadata["start_position"] >= int(start))
            & (gene_metadata["end_position"] <= int(end))
        ]["hgnc_symbol"].drop_duplicates()
        top_expr = _union_preserve(top_expr, [gene for gene in region if gene in common_genes])

    meta = gene_metadata[gene_metadata["hgnc_symbol"].isin(top_expr)].copy()
    meta["chr_arm_full"] = meta["chromosome_num"].astype(str) + meta["chr_arm"].astype(str)
    genes_by_arm = {name: list(group["hgnc_symbol"]) for name, group in meta.groupby("chr_arm_full", sort=True)}

    # Preserve fastCNV's current behavior: this loop is over bare arms p/q, not chromosome-specific arms.
    for arm in pd.unique(gene_metadata["chr_arm"]):
        if arm not in genes_by_arm:
            genes_by_arm[str(arm)] = []
        if len(genes_by_arm[str(arm)]) < 200:
            remaining = [gene for gene in common_genes if gene not in set(genes_by_arm[str(arm)])]
            top_arm = list(average_expression.reindex(remaining).sort_values(ascending=False, kind="mergesort").index[:200])
            genes_by_arm[str(arm)] = _unique_preserve([*genes_by_arm[str(arm)], *top_arm])

    if chr_arms_to_force is not None:
        arms = [chr_arms_to_force] if isinstance(chr_arms_to_force, str) else list(chr_arms_to_force)
        for arm in arms:
            forced = gene_metadata[
                gene_metadata["chromosome_num"].astype(str) + gene_metadata["chr_arm"].astype(str) == str(arm)
            ]["hgnc_symbol"]
            genes_by_arm[str(arm)] = list(forced)

    final = []
    for key in sorted(genes_by_arm):
        final.extend(genes_by_arm[key])
    final = [gene for gene in final if gene in common_genes]

    top_meta = gene_metadata[gene_metadata["hgnc_symbol"].isin(top_expr)].copy()
    top_meta = top_meta.sort_values(["chromosome_num", "start_position"], kind="mergesort")
    return SelectedGenes(top_expr_genes=top_expr, final_selected_genes=final, top_gene_metadata=top_meta)


def _scale_factor(norm: np.ndarray, columns: pd.Index, reference_cells, scale_on_reference: bool) -> np.ndarray:
    if not scale_on_reference:
        return np.nanmean(norm, axis=1)
    col_pos = {name: i for i, name in enumerate(columns)}
    if isinstance(reference_cells, dict):
        factors = []
        for cells in reference_cells.values():
            idx = [col_pos[cell] for cell in cells if cell in col_pos]
            if idx:
                factors.append(np.nanmean(norm[:, idx], axis=1))
        if not factors:
            return np.nanmean(norm, axis=1)
        return np.nanmedian(np.vstack(factors), axis=0)
    idx = [col_pos[cell] for cell in reference_cells if cell in col_pos]
    if not idx:
        return np.nanmean(norm, axis=1)
    return np.nanmean(norm[:, idx], axis=1)


def _build_genomic_windows(gene_metadata: pd.DataFrame, *, window_size: int, window_step: int) -> dict[str, list[str]]:
    windows: dict[str, list[str]] = {}
    iter_half = int(np.round(window_size / 2.0))
    for chrom in range(1, 24):
        genes_c = gene_metadata[gene_metadata["chromosome_num"] == chrom]
        for arm in pd.unique(genes_c["chr_arm"]):
            genes_arm = genes_c[genes_c["chr_arm"] == arm]
            genes = list(genes_arm["hgnc_symbol"].astype(str))
            n = len(genes)
            if n == 0:
                continue
            if n > window_size:
                starts = range(iter_half, n - iter_half, window_step)
                for idx, center in enumerate(starts, start=1):
                    name = f"{chrom}.{arm}{idx}"
                    windows[name] = genes[center - iter_half : center + iter_half + 1]
            else:
                windows[f"{chrom}.{arm}1"] = genes
    return windows


def _compute_genomic_scores(norm_counts: pd.DataFrame, windows: dict[str, list[str]]) -> np.ndarray:
    return _compute_genomic_scores_from_array(norm_counts.to_numpy(dtype=np.float64, copy=False), norm_counts.index, windows)


def _compute_genomic_scores_from_array(
    values: np.ndarray, gene_index: pd.Index, windows: dict[str, list[str]]
) -> np.ndarray:
    scores = np.empty((values.shape[1], len(windows)), dtype=np.float64)
    first_pos: dict[str, int] = {}
    for pos, gene in enumerate(gene_index.astype(str)):
        first_pos.setdefault(gene, pos)
    for j, genes in enumerate(windows.values()):
        positions = [first_pos[gene] for gene in genes if gene in first_pos]
        if len(positions) == 1:
            scores[:, j] = values[positions[0], :]
        else:
            scores[:, j] = values[positions, :].mean(axis=0)
    return scores


def _threshold_scores_like_fastcnv(
    *,
    raw_scores: np.ndarray,
    cell_names: pd.Index,
    obs: pd.DataFrame | None,
    reference_var: str | None,
    reference_cells,
    scale_on_reference: bool,
    threshold_percentile: float,
) -> np.ndarray:
    if scale_on_reference:
        if isinstance(reference_cells, dict):
            ref_cells = [cell for values in reference_cells.values() for cell in values]
        else:
            ref_cells = list(reference_cells)
        pos = {cell: idx for idx, cell in enumerate(cell_names)}
        idx = [pos[cell] for cell in ref_cells if cell in pos]
        if idx:
            return _threshold_scores(raw_scores, threshold_percentile, reference_rows=idx)
    if reference_var is not None and obs is not None and reference_var in obs.columns:
        low_values = []
        high_values = []
        aligned = obs.reindex(cell_names)
        for _, group in aligned.groupby(reference_var, sort=True):
            idx = [cell_names.get_loc(cell) for cell in group.index if cell in cell_names]
            if idx:
                q = np.nanquantile(raw_scores[idx, :], [threshold_percentile, 1.0 - threshold_percentile], axis=0)
                low_values.append(q[0])
                high_values.append(q[1])
        if low_values and high_values:
            low = np.nanmedian(np.vstack(low_values), axis=0)
            high = np.nanmedian(np.vstack(high_values), axis=0)
            trimmed = raw_scores.copy()
            inside = (trimmed > low[None, :]) & (trimmed < high[None, :])
            trimmed[inside] = 0.0
            return trimmed.T
    return _threshold_scores(raw_scores, threshold_percentile)


def _threshold_scores(
    raw_scores: np.ndarray, threshold_percentile: float, *, reference_rows: list[int] | None = None
) -> np.ndarray:
    source = raw_scores if reference_rows is None else raw_scores[reference_rows, :]
    q = np.nanquantile(source, [threshold_percentile, 1.0 - threshold_percentile], axis=0)
    trimmed = raw_scores.copy()
    inside = (trimmed >= q[0][None, :]) & (trimmed <= q[1][None, :])
    trimmed[inside] = 0.0
    return trimmed.T


def cnv_per_chromosome_arm(genomic_scores: pd.DataFrame) -> pd.DataFrame:
    arms = chromosome_arm_order()
    info = pd.Series([window_chrom_arm(name) for name in genomic_scores.index], index=genomic_scores.index)
    rows = {}
    for arm in arms:
        windows = list(info.index[info == arm])
        if windows:
            rows[f"{arm}_CNV"] = genomic_scores.loc[windows].mean(axis=0)
        else:
            rows[f"{arm}_CNV"] = pd.Series(0.0, index=genomic_scores.columns)
    return pd.DataFrame(rows).T


def classify_arm_cnv(arm_cnv: pd.DataFrame, peaks: tuple[float, float, float]) -> pd.DataFrame:
    values = arm_cnv.to_numpy()
    classes = np.full(values.shape, "no_alteration", dtype=object)
    classes[values < peaks[0]] = "loss"
    classes[values > peaks[2]] = "gain"
    return pd.DataFrame(classes, index=[f"{idx}_classification" for idx in arm_cnv.index], columns=arm_cnv.columns)


def cluster_cnv(genomic_scores: pd.DataFrame, *, k: int | None = None, h: float | None = None) -> pd.Series:
    mat = genomic_scores.T.to_numpy(dtype=np.float64)
    cells = genomic_scores.columns
    if mat.shape[0] < 2:
        return pd.Series(1, index=cells)
    dist = pdist(mat, metric="cityblock")
    hc = linkage(dist, method="ward")
    if k is None:
        k = _elbow_k(hc, dist, n_obs=mat.shape[0])
    if h is not None and k is None:
        labels = fcluster(hc, t=h, criterion="distance")
    else:
        labels = fcluster(hc, t=int(k), criterion="maxclust")
    return pd.Series(labels.astype(int), index=cells)


def merge_cnv_clusters(cell_meta: pd.DataFrame, *, merge_threshold: float = 0.98) -> pd.Series:
    if "cnv_clusters" not in cell_meta.columns:
        return pd.Series(pd.NA, index=cell_meta.index, dtype="Int64")
    arm_cols = [col for col in cell_meta.columns if col.endswith("_CNV")]
    clusters = cell_meta["cnv_clusters"].dropna().astype(int)
    if clusters.nunique() < 2 or not arm_cols:
        return cell_meta["cnv_clusters"]
    profiles = cell_meta.loc[clusters.index, arm_cols].groupby(clusters).mean()
    profiles = profiles.loc[:, (profiles != 0).any(axis=0)]
    if profiles.shape[0] < 2 or profiles.shape[1] == 0:
        return cell_meta["cnv_clusters"]
    corr = np.corrcoef(profiles.to_numpy())
    names = list(profiles.index.astype(int))
    groups = {name: name for name in names}
    for i, left in enumerate(names[:-1]):
        for j, right in enumerate(names[i + 1 :], start=i + 1):
            if corr[i, j] > merge_threshold:
                gmin, gmax = sorted((groups[left], groups[right]))
                for key, value in list(groups.items()):
                    if value == gmax:
                        groups[key] = gmin
    mapped = cell_meta["cnv_clusters"].map(lambda x: groups.get(int(x), int(x)) if pd.notna(x) else pd.NA)
    unique = {value: i + 1 for i, value in enumerate(sorted({int(v) for v in mapped.dropna()}))}
    return mapped.map(lambda x: unique[int(x)] if pd.notna(x) else pd.NA).astype("Int64")


def _elbow_k(hc, dist, *, n_obs: int) -> int:
    max_k = min(15, n_obs)
    dist_full = squareform(dist)
    k_values = np.arange(1, max_k + 1)
    wss = []
    for k in k_values:
        labels = fcluster(hc, t=int(k), criterion="maxclust")
        total = 0.0
        for cl in np.unique(labels):
            idx = np.flatnonzero(labels == cl)
            if len(idx) > 0:
                total += np.sum(dist_full[np.ix_(idx, idx)] ** 2) / (2.0 * len(idx))
        wss.append(total)
    wss = np.asarray(wss)
    if len(k_values) < 3 or np.nanmax(wss) == np.nanmin(wss):
        return int(k_values[max(0, len(k_values) // 2)])
    x_norm = (k_values - k_values.min()) / (k_values.max() - k_values.min())
    y_norm = (wss - wss.min()) / (wss.max() - wss.min())
    slopes = np.diff(y_norm) / np.diff(x_norm)
    changes = np.diff(slopes)
    if changes.size == 0:
        return 1
    threshold = 0.05 * np.nanmax(np.abs(changes))
    significant = np.flatnonzero(np.abs(changes) > threshold)
    if significant.size == 0:
        return int(k_values[int(np.ceil(len(k_values) / 2)) - 1])
    first_quarter = int(np.ceil(len(k_values) / 4))
    after = significant[significant + 1 > first_quarter]
    point = after[0] if after.size else significant[-1]
    return int(k_values[point + 1])


def _window_metadata(windows: dict[str, list[str]]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "window": list(windows.keys()),
            "chrom_arm": [window_chrom_arm(name) for name in windows.keys()],
            "n_genes": [len(genes) for genes in windows.values()],
            "genes": [";".join(genes) for genes in windows.values()],
        }
    )


def _union_preserve(left: list[str], right: list[str]) -> list[str]:
    return _unique_preserve([*left, *right])


def _unique_preserve(values: list[str]) -> list[str]:
    seen = set()
    out = []
    for value in values:
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out
