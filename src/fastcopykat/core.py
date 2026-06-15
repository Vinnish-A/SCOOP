from __future__ import annotations

from dataclasses import asdict, dataclass
from time import perf_counter

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.ndimage import uniform_filter1d
from scipy.spatial.distance import pdist

from .config import FastCopyKatConfig
from .reference import build_bins_from_annotation, filter_copykat_chromosomes, normalize_bins, normalize_gene_annotation


@dataclass(frozen=True)
class FastCopyKatResult:
    prediction: pd.DataFrame
    cna: pd.DataFrame
    gene_cna: pd.DataFrame
    cell_scores: pd.DataFrame
    manifest: dict


def run_fastcopykat(
    counts: pd.DataFrame,
    gene_annotation: pd.DataFrame,
    *,
    bins: pd.DataFrame | None = None,
    normal_cell_names: list[str] | tuple[str, ...] | None = None,
    sample_name: str = "sample",
    config: FastCopyKatConfig | None = None,
) -> FastCopyKatResult:
    """Run an independent Python approximation of CopyKAT's main pipeline.

    The input matrix follows CopyKAT's convention: genes are rows and cells are columns.
    """

    cfg = config or FastCopyKatConfig()
    cfg.validate()
    timings: dict[str, float] = {}
    total_start = perf_counter()

    start = perf_counter()
    annotation = normalize_gene_annotation(gene_annotation)
    annotation = filter_copykat_chromosomes(annotation, genome=cfg.genome)
    bin_index = normalize_bins(bins) if bins is not None else build_bins_from_annotation(annotation, cfg.bin_size)
    bin_index = filter_copykat_chromosomes(bin_index, genome=cfg.genome)
    matrix, annotation = _align_and_filter_counts(counts, annotation, cfg)
    timings["filter_seconds"] = perf_counter() - start

    start = perf_counter()
    count_values = matrix.to_numpy(dtype=np.float64, copy=False)
    transformed = np.log(np.sqrt(count_values) + np.sqrt(count_values + 1.0))
    transformed -= np.nanmean(transformed, axis=0, keepdims=True)
    smoothed = _smooth_by_chromosome(transformed, annotation["chrom"].to_numpy(), cfg.window_size)
    timings["smooth_seconds"] = perf_counter() - start

    start = perf_counter()
    baseline, baseline_cells, baseline_labels = _estimate_baseline(
        smoothed,
        cell_names=np.asarray(matrix.columns, dtype=object),
        normal_cell_names=normal_cell_names or (),
        cfg=cfg,
    )
    relative = smoothed - baseline[:, None]
    upper_dr = cfg.upper_detection_rate
    if relative.shape[0] < 7000:
        upper_dr = cfg.low_detection_rate
    keep_upper = (count_values > 0).mean(axis=1) >= upper_dr
    relative = relative[keep_upper]
    seg_annotation = annotation.loc[keep_upper].reset_index(drop=True)
    timings["baseline_seconds"] = perf_counter() - start

    start = perf_counter()
    gene_cna_values = _segment_profiles(relative, seg_annotation["chrom"].to_numpy(), cfg)
    timings["segment_seconds"] = perf_counter() - start

    start = perf_counter()
    bin_values = _genes_to_bins(gene_cna_values, seg_annotation, bin_index)
    adjusted_bins, prediction, cell_scores = _predict_cells(
        bin_values,
        bin_chrom=bin_index["chrom"].to_numpy(),
        cell_names=np.asarray(matrix.columns, dtype=object),
        baseline_cells=baseline_cells,
        cfg=cfg,
    )
    timings["bin_and_predict_seconds"] = perf_counter() - start

    gene_cna = pd.DataFrame(gene_cna_values, columns=matrix.columns)
    gene_cna.insert(0, "abspos", seg_annotation["abspos"].to_numpy())
    gene_cna.insert(0, "chrompos", seg_annotation["chrompos"].to_numpy())
    gene_cna.insert(0, "chrom", seg_annotation["chrom"].to_numpy())

    cna = pd.DataFrame(adjusted_bins, columns=matrix.columns)
    cna.insert(0, "abspos", bin_index["abspos"].to_numpy())
    cna.insert(0, "chrompos", bin_index["chrompos"].to_numpy())
    cna.insert(0, "chrom", bin_index["chrom"].to_numpy())

    timings["total_seconds"] = perf_counter() - total_start
    manifest = {
        "schema_version": "fastcopykat.manifest.v1",
        "sample_name": sample_name,
        "engine": "fastcopykat",
        "config": asdict(cfg),
        "n_input_genes": int(counts.shape[0]),
        "n_input_cells": int(counts.shape[1]),
        "n_retained_genes": int(gene_cna_values.shape[0]),
        "n_retained_cells": int(matrix.shape[1]),
        "n_bins": int(bin_values.shape[0]),
        "n_baseline_cells": int(len(baseline_cells)),
        "baseline_cells": list(map(str, baseline_cells[:100])),
        "baseline_cluster_labels": baseline_labels,
        "prediction_summary": {
            "counts": prediction["copykat.pred"].value_counts().to_dict(),
            "thresholds": {
                key: float(cell_scores[key].iloc[0])
                for key in ("cnv_score_threshold", "chrom_burden_threshold")
                if key in cell_scores.columns and len(cell_scores) > 0
            },
        },
        "timings": {key: round(value, 6) for key, value in timings.items()},
    }
    return FastCopyKatResult(
        prediction=prediction,
        cna=cna,
        gene_cna=gene_cna,
        cell_scores=cell_scores,
        manifest=manifest,
    )


def _align_and_filter_counts(
    counts: pd.DataFrame, annotation: pd.DataFrame, cfg: FastCopyKatConfig
) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = counts.copy()
    counts.index = counts.index.astype(str)
    counts.columns = counts.columns.astype(str)
    annotation = annotation[annotation["gene"].isin(counts.index)].copy()
    if annotation.empty:
        raise ValueError("no genes overlap between counts and gene annotation")
    counts = counts.loc[annotation["gene"]]
    detected_per_cell = (counts > 0).sum(axis=0)
    cell_keep = detected_per_cell >= cfg.min_gene_per_cell
    counts = counts.loc[:, cell_keep]
    if counts.shape[1] < 2:
        raise ValueError("fewer than two cells remain after min_gene_per_cell filtering")
    detection_rate = (counts > 0).mean(axis=1)
    gene_keep = detection_rate >= cfg.low_detection_rate
    counts = counts.loc[gene_keep]
    annotation = annotation.loc[gene_keep.to_numpy()].reset_index(drop=True)

    chrom_counts = (counts > 0).groupby(annotation["chrom"].to_numpy()).sum()
    valid_cells = (chrom_counts >= cfg.min_gene_per_chromosome).all(axis=0)
    counts = counts.loc[:, valid_cells]
    if counts.shape[1] < 2:
        raise ValueError("fewer than two cells remain after chromosome coverage filtering")
    return counts.astype("float64", copy=False), annotation


def _smooth_by_chromosome(values: np.ndarray, chrom: np.ndarray, window_size: int) -> np.ndarray:
    smoothed = np.empty_like(values, dtype=np.float64)
    for _, idx in _chromosome_slices(chrom):
        segment = values[idx, :]
        width = min(window_size, max(3, segment.shape[0]))
        if width % 2 == 0:
            width -= 1
        smoothed[idx, :] = uniform_filter1d(segment, size=width, axis=0, mode="nearest")
    return smoothed


def _estimate_baseline(
    values: np.ndarray,
    *,
    cell_names: np.ndarray,
    normal_cell_names: list[str] | tuple[str, ...],
    cfg: FastCopyKatConfig,
) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    normal_set = set(map(str, normal_cell_names))
    normal_idx = np.asarray([name in normal_set for name in cell_names])
    if normal_idx.any():
        return np.nanmedian(values[:, normal_idx], axis=1), cell_names[normal_idx], {"known_normals": int(normal_idx.sum())}

    n_cells = values.shape[1]
    max_k = min(cfg.max_baseline_clusters, max(2, n_cells // cfg.min_cluster_cells))
    if n_cells <= max(3, cfg.min_cluster_cells * 2):
        score = np.nanmedian(np.abs(values), axis=0)
        keep = score <= np.quantile(score, cfg.prediction_score_quantile)
        return np.nanmedian(values[:, keep], axis=1), cell_names[keep], {"small_dataset": int(keep.sum())}

    distance = pdist(values.T, metric=cfg.distance)
    method = "ward" if cfg.distance == "euclidean" else "average"
    tree = linkage(distance, method=method)
    best_labels = np.ones(n_cells, dtype=int)
    best_k = 1
    for k in range(max_k, 1, -1):
        labels = fcluster(tree, k, criterion="maxclust")
        sizes = pd.Series(labels).value_counts()
        if sizes.min() >= cfg.min_cluster_cells:
            best_labels = labels
            best_k = k
            break
    scores = {}
    for label in np.unique(best_labels):
        idx = best_labels == label
        med = np.nanmedian(values[:, idx], axis=1)
        scores[int(label)] = float(np.nanmedian(np.abs(med)) + np.nanstd(med))
    baseline_label = min(scores, key=scores.get)
    keep = best_labels == baseline_label
    return (
        np.nanmedian(values[:, keep], axis=1),
        cell_names[keep],
        {"k": int(best_k), "baseline_label": int(baseline_label), "baseline_cells": int(keep.sum())},
    )


def _segment_profiles(values: np.ndarray, chrom: np.ndarray, cfg: FastCopyKatConfig) -> np.ndarray:
    segmented = np.empty_like(values, dtype=np.float64)
    for _, idx in _chromosome_slices(chrom):
        block = values[idx, :]
        if block.shape[0] < cfg.window_size * 2:
            segmented[idx, :] = np.nanmedian(block, axis=0, keepdims=True)
            continue
        width = min(cfg.window_size, block.shape[0])
        smooth = uniform_filter1d(block, size=width, axis=0, mode="nearest")
        window = max(3, width)
        left = uniform_filter1d(smooth, size=window, axis=0, mode="nearest")
        delta = np.nanmedian(np.abs(np.diff(left, axis=0)), axis=1)
        cutoff = max(cfg.segmentation_threshold, float(np.nanquantile(delta, 0.90)))
        breaks = np.flatnonzero(delta >= cutoff) + 1
        points = np.unique(np.concatenate(([0], breaks, [block.shape[0]])))
        base = 0 if idx.start is None else idx.start
        for start, end in zip(points[:-1], points[1:]):
            if end <= start:
                continue
            segmented[base + start : base + end, :] = np.nanmedian(block[start:end, :], axis=0, keepdims=True)
    return segmented


def _genes_to_bins(values: np.ndarray, annotation: pd.DataFrame, bins: pd.DataFrame) -> np.ndarray:
    out = np.full((bins.shape[0], values.shape[1]), np.nan, dtype=np.float64)
    for chrom, bin_idx in bins.groupby("chrom", sort=False).groups.items():
        gene_idx = np.flatnonzero(annotation["chrom"].to_numpy() == str(chrom))
        if gene_idx.size == 0:
            continue
        gene_pos = annotation.iloc[gene_idx]["chrompos"].to_numpy(dtype=float)
        bin_pos = bins.iloc[list(bin_idx)]["chrompos"].to_numpy(dtype=float)
        order = np.argsort(gene_pos)
        gene_pos = gene_pos[order]
        gene_values = values[gene_idx[order], :]
        edges = np.r_[0, np.searchsorted(gene_pos, (bin_pos[:-1] + bin_pos[1:]) / 2.0), gene_pos.size]
        rows = list(bin_idx)
        for local_i, row in enumerate(rows):
            start, end = int(edges[local_i]), int(edges[local_i + 1])
            if end > start:
                out[row, :] = np.nanmedian(gene_values[start:end, :], axis=0)
    return _fill_missing_bins(out)


def _fill_missing_bins(values: np.ndarray) -> np.ndarray:
    out = values.copy()
    row_ok = ~np.isnan(out).all(axis=1)
    if not row_ok.any():
        return np.nan_to_num(out, nan=0.0)
    valid = np.flatnonzero(row_ok)
    missing = np.flatnonzero(~row_ok)
    for row in missing:
        nearest = valid[np.argmin(np.abs(valid - row))]
        out[row, :] = out[nearest, :]
    return np.nan_to_num(out, nan=0.0)


def _predict_cells(
    bins: np.ndarray,
    *,
    bin_chrom: np.ndarray,
    cell_names: np.ndarray,
    baseline_cells: np.ndarray,
    cfg: FastCopyKatConfig,
) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    baseline_set = set(map(str, baseline_cells))
    baseline_idx = np.asarray([name in baseline_set for name in cell_names])
    if not baseline_idx.any():
        score0 = np.nanmedian(np.abs(bins), axis=0)
        baseline_idx = score0 <= np.quantile(score0, cfg.prediction_score_quantile)
    baseline_profile = np.nanmedian(bins[:, baseline_idx], axis=1)
    adjusted = bins - baseline_profile[:, None]
    scores = np.nanmedian(np.abs(adjusted), axis=0)
    chrom_burden = _chromosome_burden(adjusted, bin_chrom, top_n=cfg.chromosome_rescue_top_n)
    normal_scores = scores[baseline_idx]
    normal_chrom_burden = chrom_burden[baseline_idx]
    center = float(np.nanmedian(normal_scores))
    spread = float(1.4826 * np.nanmedian(np.abs(normal_scores - center)))
    threshold = center + cfg.prediction_mad_multiplier * spread
    chrom_center = float(np.nanmedian(normal_chrom_burden))
    chrom_spread = float(1.4826 * np.nanmedian(np.abs(normal_chrom_burden - chrom_center)))
    chrom_threshold = chrom_center + cfg.chromosome_rescue_mad_multiplier * chrom_spread
    labels = np.where((scores > threshold) | (chrom_burden > chrom_threshold), "aneuploid", "diploid")
    labels[baseline_idx] = "diploid"
    prediction = pd.DataFrame({"cell.names": cell_names.astype(str), "copykat.pred": labels})
    cell_scores = pd.DataFrame(
        {
            "cell.names": cell_names.astype(str),
            "copykat.pred": labels,
            "cnv_score": scores,
            "chrom_burden_score": chrom_burden,
            "is_baseline_cell": baseline_idx,
            "cnv_score_threshold": threshold,
            "chrom_burden_threshold": chrom_threshold,
        }
    )
    return adjusted, prediction, cell_scores


def _chromosome_burden(adjusted_bins: np.ndarray, bin_chrom: np.ndarray, *, top_n: int) -> np.ndarray:
    per_chrom = []
    chrom_values = bin_chrom.astype(str)
    for chrom in pd.unique(pd.Series(chrom_values)):
        idx = chrom_values == str(chrom)
        if idx.any():
            per_chrom.append(np.nanmedian(np.abs(adjusted_bins[idx, :]), axis=0))
    if not per_chrom:
        return np.nanmedian(np.abs(adjusted_bins), axis=0)
    burden = np.vstack(per_chrom)
    n = min(top_n, burden.shape[0])
    top = np.partition(burden, kth=burden.shape[0] - n, axis=0)[-n:, :]
    return np.nanmean(top, axis=0)


def _chromosome_slices(chrom: np.ndarray):
    start = 0
    for i in range(1, len(chrom) + 1):
        if i == len(chrom) or chrom[i] != chrom[start]:
            yield chrom[start], slice(start, i)
            start = i
