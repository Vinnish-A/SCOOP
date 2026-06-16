from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from .abundance_data import AbundanceTable, build_sample_celltype_counts_from_h5ad, read_sample_celltype_counts, write_abundance_inputs
from .abundance_design import (
    build_feature_design,
    encode_binary_labels,
    encode_multiclass_labels,
    encode_survival,
    merge_metadata,
    parse_covariates,
)
from .abundance_loss import cox_partial_nll
from .abundance_metrics import (
    benjamini_hochberg,
    binary_metrics,
    concordance_index,
    correlation_z_p,
    multiclass_metrics,
    welch_z_p,
)
from .abundance_result import AbundanceResult
from .abundance_train import AbundanceTrainer, AbundanceTrainingConfig
from .abundance_mil import run_mil_abundance


def load_abundance_table(
    input_h5ad: str | Path | None = None,
    counts: str | Path | None = None,
    sample_key: str = "sample_id",
    celltype_key: str = "cell_type_lvl3",
    min_cells_per_sample: int = 20,
    min_total_cells_per_celltype: int = 50,
) -> AbundanceTable:
    if input_h5ad is None and counts is None:
        raise ValueError("provide either input_h5ad or counts")
    if input_h5ad is not None:
        adata = ad.read_h5ad(input_h5ad)
        return build_sample_celltype_counts_from_h5ad(
            adata,
            sample_key=sample_key,
            celltype_key=celltype_key,
            min_cells_per_sample=min_cells_per_sample,
            min_total_cells_per_celltype=min_total_cells_per_celltype,
        )
    table = read_sample_celltype_counts(counts)
    table.sample_key = sample_key
    table.celltype_key = celltype_key
    return table


def _read_metadata(path: str | Path | None) -> pd.DataFrame | None:
    if path is None:
        return None
    meta = pd.read_csv(path, sep="\t")
    return meta


def _subset_table(table: AbundanceTable, sample_index: pd.Index) -> AbundanceTable:
    sample_index = pd.Index(sample_index.astype(str))
    return AbundanceTable(
        counts=table.counts.loc[sample_index],
        proportions=table.proportions.loc[sample_index],
        metadata=table.metadata.loc[sample_index],
        sample_key=table.sample_key,
        celltype_key=table.celltype_key,
    )


def _training_config(kwargs: dict[str, Any]) -> AbundanceTrainingConfig:
    return AbundanceTrainingConfig(
        learning_rate=float(kwargs.get("learning_rate", 1e-3)),
        weight_decay=float(kwargs.get("weight_decay", 1e-4)),
        max_epochs=int(kwargs.get("max_epochs", 500)),
        random_seed=int(kwargs.get("random_seed", 0)),
        binary_class_weight=kwargs.get("binary_class_weight", "balanced"),
    )


def _manifest(
    *,
    mode: str,
    table: AbundanceTable,
    sample_key: str,
    celltype_key: str,
    transform: str,
    n_events: int | None = None,
    n_classes: int | None = None,
    metrics: dict[str, Any] | None = None,
    model: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "schema_version": "fastde.abundance.v1",
        "engine": "fastde.abundance",
        "scsurvival_compat": True,
        "scsurvival_source": {
            "repository": "https://github.com/cliffren/scSurvival",
            "commit": "a76de0a00035e4c4d49df4a06001b392c8014105",
            "version": "1.3.0",
            "license": "GPL-3.0",
        },
        "mode": mode,
        "sample_key": sample_key,
        "celltype_key": celltype_key,
        "n_samples": table.n_samples,
        "n_celltypes": table.n_celltypes,
        "n_events": n_events,
        "n_classes": n_classes,
        "transform": transform,
        "model": model or {},
        "loss": {"task": mode, "count_auxiliary": "documented_not_enabled_in_numpy_backend"},
        "metrics": metrics or {},
        "outputs": {},
        "diagnostics": {},
    }


def _coef_for_abundance(model, n_abundance: int, class_index: int | None = None, reference_index: int | None = None) -> np.ndarray:
    coef = np.asarray(model.coef_, dtype=float)
    if coef.ndim == 1:
        return coef[:n_abundance]
    if class_index is None:
        return coef.mean(axis=0)[:n_abundance]
    base = coef[class_index, :n_abundance]
    if reference_index is not None:
        base = base - coef[reference_index, :n_abundance]
    return base


def _binary_results(
    table: AbundanceTable,
    labels: np.ndarray,
    scores: np.ndarray,
    model,
    positive_label: str,
    negative_label: str,
    mode: str,
) -> pd.DataFrame:
    rows = []
    props = table.proportions
    coef = _coef_for_abundance(model, props.shape[1])
    for i, cell_type in enumerate(props.columns):
        pos = props.iloc[labels == 1, i].to_numpy(dtype=float)
        neg = props.iloc[labels == 0, i].to_numpy(dtype=float)
        z, pvalue = welch_z_p(pos, neg)
        c = float(coef[i])
        rows.append(
            {
                "cell_type": cell_type,
                "coef_or_importance": c,
                "direction": "positive" if c >= 0 else "negative",
                "odds_ratio_approx": float(np.exp(np.clip(c, -20, 20))),
                "mean_prop_positive": float(pos.mean()) if len(pos) else np.nan,
                "mean_prop_negative": float(neg.mean()) if len(neg) else np.nan,
                "delta_prop": float(pos.mean() - neg.mean()) if len(pos) and len(neg) else np.nan,
                "z": z,
                "pvalue": pvalue,
                "model": "scsurvival_like_abundance_logistic",
                "mode": mode,
                "positive_label": positive_label,
                "negative_label": negative_label,
            }
        )
    out = pd.DataFrame(rows)
    out["padj"] = benjamini_hochberg(out["pvalue"].to_numpy(dtype=float))
    return out.sort_values(["padj", "pvalue", "cell_type"], na_position="last").reset_index(drop=True)


def _survival_results(table: AbundanceTable, time: np.ndarray, event: np.ndarray, risk: np.ndarray, model) -> pd.DataFrame:
    rows = []
    props = table.proportions
    coef = _coef_for_abundance(model, props.shape[1])
    high = risk >= np.median(risk)
    cindex = concordance_index(time, risk, event)
    for i, cell_type in enumerate(props.columns):
        feature = props.iloc[:, i].to_numpy(dtype=float)
        z, pvalue = correlation_z_p(feature, risk)
        c = float(coef[i])
        rows.append(
            {
                "cell_type": cell_type,
                "coef_or_importance": c,
                "hazard_ratio_approx": float(np.exp(np.clip(c, -20, 20))),
                "direction": "high_risk" if c >= 0 else "low_risk",
                "mean_prop_high_risk": float(feature[high].mean()) if high.any() else np.nan,
                "mean_prop_low_risk": float(feature[~high].mean()) if (~high).any() else np.nan,
                "z": z,
                "pvalue": pvalue,
                "concordance_index": cindex,
                "n_events": int(event.sum()),
                "n_samples": int(len(time)),
                "model": "scsurvival_like_abundance_cox",
                "mode": "survival",
            }
        )
    out = pd.DataFrame(rows)
    out["padj"] = benjamini_hochberg(out["pvalue"].to_numpy(dtype=float))
    return out.sort_values(["padj", "pvalue", "cell_type"], na_position="last").reset_index(drop=True)


def _multiclass_results(table: AbundanceTable, labels: np.ndarray, classes: list[str], reference_idx: int, model) -> pd.DataFrame:
    rows = []
    props = table.proportions
    for class_idx, cls in enumerate(classes):
        if class_idx == reference_idx:
            continue
        coef = _coef_for_abundance(model, props.shape[1], class_index=class_idx, reference_index=reference_idx)
        for i, cell_type in enumerate(props.columns):
            cur = props.iloc[labels == class_idx, i].to_numpy(dtype=float)
            ref = props.iloc[labels == reference_idx, i].to_numpy(dtype=float)
            z, pvalue = welch_z_p(cur, ref)
            c = float(coef[i])
            rows.append(
                {
                    "cell_type": cell_type,
                    "class_or_contrast": f"{cls}_vs_{classes[reference_idx]}",
                    "coef_or_importance": c,
                    "direction": cls if c >= 0 else classes[reference_idx],
                    "mean_prop_class": float(cur.mean()) if len(cur) else np.nan,
                    "mean_prop_reference": float(ref.mean()) if len(ref) else np.nan,
                    "delta_prop": float(cur.mean() - ref.mean()) if len(cur) and len(ref) else np.nan,
                    "z": z,
                    "pvalue": pvalue,
                    "model": "scsurvival_like_abundance_multinomial",
                    "mode": "multiclass",
                }
            )
    out = pd.DataFrame(rows)
    out["padj"] = benjamini_hochberg(out["pvalue"].to_numpy(dtype=float))
    return out.sort_values(["padj", "pvalue", "cell_type"], na_position="last").reset_index(drop=True)


def _continuous_results(table: AbundanceTable, values: np.ndarray, scores: np.ndarray, model) -> pd.DataFrame:
    rows = []
    props = table.proportions
    coef = _coef_for_abundance(model, props.shape[1])
    for i, cell_type in enumerate(props.columns):
        z, pvalue = correlation_z_p(props.iloc[:, i].to_numpy(dtype=float), values)
        c = float(coef[i])
        rows.append(
            {
                "cell_type": cell_type,
                "coef_or_importance": c,
                "direction": "positive" if c >= 0 else "negative",
                "z": z,
                "pvalue": pvalue,
                "padj": np.nan,
                "model": "scsurvival_like_abundance_ridge",
                "mode": "continuous",
            }
        )
    out = pd.DataFrame(rows)
    out["padj"] = benjamini_hochberg(out["pvalue"].to_numpy(dtype=float))
    return out.sort_values(["padj", "pvalue", "cell_type"], na_position="last").reset_index(drop=True)


def run_abundance(
    *,
    mode: str,
    output_dir: str | Path,
    input_h5ad: str | Path | None = None,
    counts: str | Path | None = None,
    metadata: str | Path | pd.DataFrame | None = None,
    sample_key: str = "sample_id",
    celltype_key: str = "cell_type_lvl3",
    label_col: str | None = None,
    positive_label: str | None = None,
    negative_label: str | None = None,
    reference_level: str | None = None,
    time_col: str | None = None,
    event_col: str | None = None,
    value_col: str | None = None,
    covariates: str | list[str] | None = None,
    transform: str = "clr",
    pseudocount: float = 0.5,
    min_cells_per_sample: int = 20,
    min_total_cells_per_celltype: int = 50,
    **kwargs: Any,
) -> AbundanceResult:
    mode = mode.lower()
    if mode == "condition":
        mode = "binary"
    model_backend = str(kwargs.get("model_backend", "scsurvival_mil"))
    table = load_abundance_table(
        input_h5ad=input_h5ad,
        counts=counts,
        sample_key=sample_key,
        celltype_key=celltype_key,
        min_cells_per_sample=min_cells_per_sample,
        min_total_cells_per_celltype=min_total_cells_per_celltype,
    )
    if isinstance(metadata, pd.DataFrame):
        meta = metadata
    else:
        meta = _read_metadata(metadata)
    table = merge_metadata(table, meta, sample_key=sample_key)
    outdir = Path(output_dir)
    write_abundance_inputs(table, outdir)
    cov = parse_covariates(covariates)

    if model_backend in {"scsurvival_mil", "mil"}:
        return run_mil_abundance(
            mode=mode,
            table=table,
            output_dir=outdir,
            input_h5ad=input_h5ad,
            sample_key=sample_key,
            celltype_key=celltype_key,
            label_col=label_col,
            positive_label=positive_label,
            negative_label=negative_label,
            reference_level=reference_level,
            time_col=time_col,
            event_col=event_col,
            value_col=value_col,
            covariates=covariates,
            max_instances_per_sample=int(kwargs.get("max_instances_per_sample", 2000)),
            hidden_dim=int(kwargs.get("hidden_dim", 64)),
            dropout=float(kwargs.get("dropout", 0.1)),
            learning_rate=float(kwargs.get("learning_rate", 1e-3)),
            weight_decay=float(kwargs.get("weight_decay", 1e-4)),
            max_epochs=int(kwargs.get("max_epochs", 500)),
            random_seed=int(kwargs.get("random_seed", 0)),
            survival_loss=str(kwargs.get("survival_loss", "cox")),
            manifest_factory=_manifest,
        )

    if mode == "binary":
        if not label_col or positive_label is None or negative_label is None:
            raise ValueError("binary mode requires label_col, positive_label, and negative_label")
        y, keep = encode_binary_labels(table.metadata, label_col, positive_label, negative_label)
        table = _subset_table(table, keep)
        design = build_feature_design(table, transform=transform, pseudocount=pseudocount, covariates=cov)
        trainer = AbundanceTrainer("binary", _training_config(kwargs))
        model = trainer.fit(design.features, y=y)
        logit = model.predict_score(design.features.to_numpy(dtype=float))
        metrics = binary_metrics(y, logit)
        results = _binary_results(table, y.astype(int), logit, model, str(positive_label), str(negative_label), "binary")
        predictions = pd.DataFrame({"sample_id": table.counts.index, "label": y.astype(int), "logit": logit, "probability": 1 / (1 + np.exp(-logit))})
        manifest = _manifest(mode="binary", table=table, sample_key=sample_key, celltype_key=celltype_key, transform=transform, n_classes=2, metrics=metrics, model={"backend": "numpy_sklearn"})
    elif mode == "survival":
        if not time_col or not event_col:
            raise ValueError("survival mode requires time_col and event_col")
        time, event, keep = encode_survival(table.metadata, time_col, event_col)
        table = _subset_table(table, keep)
        design = build_feature_design(table, transform=transform, pseudocount=pseudocount, covariates=cov)
        trainer = AbundanceTrainer("survival", _training_config(kwargs))
        model = trainer.fit(design.features, time=time, event=event)
        risk = model.predict_score(design.features.to_numpy(dtype=float))
        metrics = {"cox_partial_nll": cox_partial_nll(risk, time, event), "concordance_index": concordance_index(time, risk, event)}
        results = _survival_results(table, time, event, risk, model)
        predictions = pd.DataFrame({"sample_id": table.counts.index, "time": time, "event": event.astype(int), "risk": risk})
        manifest = _manifest(mode="survival", table=table, sample_key=sample_key, celltype_key=celltype_key, transform=transform, n_events=int(event.sum()), metrics=metrics, model={"backend": "numpy_cox"})
    elif mode == "multiclass":
        if not label_col:
            raise ValueError("multiclass mode requires label_col")
        y, classes, keep, reference_idx = encode_multiclass_labels(table.metadata, label_col, reference_level)
        if reference_idx is None:
            reference_idx = 0
        table = _subset_table(table, keep)
        design = build_feature_design(table, transform=transform, pseudocount=pseudocount, covariates=cov)
        trainer = AbundanceTrainer("multiclass", _training_config(kwargs))
        model = trainer.fit(design.features, y=y)
        logits = model.predict_logits(design.features.to_numpy(dtype=float))
        metrics = multiclass_metrics(y, logits, classes)
        results = _multiclass_results(table, y, classes, reference_idx, model)
        predictions = pd.DataFrame({"sample_id": table.counts.index, "label": [classes[int(i)] for i in y], "predicted_label": [classes[int(i)] for i in logits.argmax(axis=1)]})
        for i, cls in enumerate(classes):
            predictions[f"prob_{cls}"] = np.exp(logits[:, i]) / np.exp(logits).sum(axis=1)
        manifest = _manifest(mode="multiclass", table=table, sample_key=sample_key, celltype_key=celltype_key, transform=transform, n_classes=len(classes), metrics=metrics, model={"backend": "numpy_sklearn"})
    elif mode == "continuous":
        if not value_col or value_col not in table.metadata.columns:
            raise ValueError("continuous mode requires value_col in metadata")
        values = pd.to_numeric(table.metadata[value_col], errors="coerce")
        keep = table.metadata.index[values.notna()]
        y = values.loc[keep].to_numpy(dtype=float)
        table = _subset_table(table, keep)
        design = build_feature_design(table, transform=transform, pseudocount=pseudocount, covariates=cov)
        trainer = AbundanceTrainer("continuous", _training_config(kwargs))
        model = trainer.fit(design.features, y=y)
        pred = model.predict_score(design.features.to_numpy(dtype=float))
        z, pvalue = correlation_z_p(pred, y)
        metrics = {"prediction_correlation_z": z, "prediction_correlation_pvalue": pvalue}
        results = _continuous_results(table, y, pred, model)
        predictions = pd.DataFrame({"sample_id": table.counts.index, "value": y, "prediction": pred})
        manifest = _manifest(mode="continuous", table=table, sample_key=sample_key, celltype_key=celltype_key, transform=transform, metrics=metrics, model={"backend": "numpy_sklearn"})
    else:
        raise ValueError(f"unknown abundance mode: {mode}")

    result = AbundanceResult(mode=mode, results=results, predictions=predictions, metrics=metrics, manifest=manifest)
    result.write(outdir)
    return result
