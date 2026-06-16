from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd

from .abundance_data import AbundanceTable


@dataclass
class FeatureDesign:
    abundance: pd.DataFrame
    covariates: pd.DataFrame
    features: pd.DataFrame
    abundance_columns: list[str]
    covariate_columns: list[str]


def parse_covariates(covariates: str | Iterable[str] | None) -> list[str]:
    if covariates is None:
        return []
    if isinstance(covariates, str):
        return [item.strip() for item in covariates.split(",") if item.strip()]
    return [str(item) for item in covariates]


def merge_metadata(table: AbundanceTable, metadata: pd.DataFrame | None, sample_key: str | None = None) -> AbundanceTable:
    if metadata is None:
        meta = table.metadata.reindex(table.counts.index).copy()
    else:
        meta = metadata.copy()
        key = sample_key or table.sample_key
        if key in meta.columns:
            meta.index = meta[key].astype(str)
        meta.index = meta.index.astype(str)
        meta = meta.reindex(table.counts.index)
    meta.index.name = table.sample_key
    return AbundanceTable(table.counts, table.proportions, meta, table.sample_key, table.celltype_key)


def abundance_transform(
    counts: pd.DataFrame,
    transform: str = "clr",
    pseudocount: float = 0.5,
) -> pd.DataFrame:
    values = counts.to_numpy(dtype=float)
    if np.any(values < 0):
        raise ValueError("abundance counts must be non-negative")
    transform = transform.lower()
    adjusted = values + float(pseudocount)
    props = adjusted / adjusted.sum(axis=1, keepdims=True)
    if transform == "proportion":
        z = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0).to_numpy(dtype=float)
    elif transform == "clr":
        logp = np.log(props)
        z = logp - logp.mean(axis=1, keepdims=True)
    elif transform == "logit_proportion":
        eps = 1e-6
        p = np.clip(props, eps, 1.0 - eps)
        z = np.log(p / (1.0 - p))
    elif transform == "arcsin_sqrt":
        z = np.arcsin(np.sqrt(props))
    elif transform == "raw_count_with_offset":
        z = np.log1p(adjusted)
    else:
        raise ValueError(f"unknown abundance transform: {transform}")
    return pd.DataFrame(z, index=counts.index, columns=counts.columns)


def build_covariate_matrix(metadata: pd.DataFrame, covariates: list[str]) -> pd.DataFrame:
    if not covariates:
        return pd.DataFrame(index=metadata.index)
    missing = [col for col in covariates if col not in metadata.columns]
    if missing:
        raise KeyError(f"metadata is missing covariates: {missing}")
    parts: list[pd.DataFrame] = []
    for col in covariates:
        series = metadata[col]
        if pd.api.types.is_numeric_dtype(series):
            values = series.astype(float)
            std = values.std(ddof=0)
            if np.isfinite(std) and std > 0:
                values = (values - values.mean()) / std
            else:
                values = values * 0.0
            parts.append(pd.DataFrame({col: values}, index=metadata.index))
        else:
            dummies = pd.get_dummies(series.astype(str), prefix=col, drop_first=True, dtype=float)
            dummies.index = metadata.index
            parts.append(dummies)
    if not parts:
        return pd.DataFrame(index=metadata.index)
    return pd.concat(parts, axis=1).astype(float)


def build_feature_design(
    table: AbundanceTable,
    transform: str = "clr",
    pseudocount: float = 0.5,
    covariates: list[str] | None = None,
) -> FeatureDesign:
    abundance = abundance_transform(table.counts, transform=transform, pseudocount=pseudocount)
    cov = build_covariate_matrix(table.metadata, covariates or [])
    features = pd.concat([abundance, cov], axis=1).astype(float)
    return FeatureDesign(
        abundance=abundance,
        covariates=cov,
        features=features,
        abundance_columns=list(abundance.columns),
        covariate_columns=list(cov.columns),
    )


def encode_binary_labels(
    metadata: pd.DataFrame,
    label_col: str,
    positive_label: str,
    negative_label: str,
) -> tuple[np.ndarray, pd.Index]:
    if label_col not in metadata.columns:
        raise KeyError(f"metadata is missing label column {label_col!r}")
    labels = metadata[label_col].astype(str)
    keep = labels.isin([str(positive_label), str(negative_label)])
    y = (labels.loc[keep] == str(positive_label)).astype(float).to_numpy()
    if len(np.unique(y)) != 2:
        raise ValueError("binary mode requires both positive and negative samples")
    return y, metadata.index[keep]


def encode_multiclass_labels(
    metadata: pd.DataFrame,
    label_col: str,
    reference_level: str | None = None,
) -> tuple[np.ndarray, list[str], pd.Index, int | None]:
    if label_col not in metadata.columns:
        raise KeyError(f"metadata is missing label column {label_col!r}")
    labels = metadata[label_col].astype(str)
    keep = labels.notna()
    classes = sorted(labels.loc[keep].unique().tolist())
    if reference_level is not None:
        ref = str(reference_level)
        if ref not in classes:
            raise ValueError(f"reference_level {ref!r} is not present")
        classes = [ref] + [cls for cls in classes if cls != ref]
    if len(classes) < 2:
        raise ValueError("multiclass mode requires at least two classes")
    mapping = {cls: i for i, cls in enumerate(classes)}
    y = labels.loc[keep].map(mapping).astype(int).to_numpy()
    ref_idx = 0 if reference_level is not None else None
    return y, classes, metadata.index[keep], ref_idx


def encode_survival(metadata: pd.DataFrame, time_col: str, event_col: str) -> tuple[np.ndarray, np.ndarray, pd.Index]:
    for col in [time_col, event_col]:
        if col not in metadata.columns:
            raise KeyError(f"metadata is missing survival column {col!r}")
    time = pd.to_numeric(metadata[time_col], errors="coerce")
    event = pd.to_numeric(metadata[event_col], errors="coerce")
    keep = time.notna() & event.notna()
    event_values = event.loc[keep].astype(int).to_numpy()
    if event_values.sum() < 1:
        raise ValueError("survival mode requires at least one observed event")
    return time.loc[keep].astype(float).to_numpy(), event_values.astype(float), metadata.index[keep]
