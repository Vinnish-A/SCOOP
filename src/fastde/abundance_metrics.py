from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import softmax


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    restored = np.empty_like(adj)
    restored[order] = np.clip(adj, 0, 1)
    out[valid] = restored
    return out


def concordance_index(time: np.ndarray, risk: np.ndarray, event: np.ndarray) -> float:
    time = np.asarray(time, dtype=float)
    risk = np.asarray(risk, dtype=float)
    event = np.asarray(event, dtype=float)
    comparable = 0
    concordant = 0.0
    for i in range(len(time)):
        if event[i] <= 0:
            continue
        for j in range(len(time)):
            if time[i] >= time[j]:
                continue
            comparable += 1
            if risk[i] > risk[j]:
                concordant += 1.0
            elif risk[i] == risk[j]:
                concordant += 0.5
    return float(concordant / comparable) if comparable else float("nan")


def binary_metrics(y: np.ndarray, logit: np.ndarray) -> dict[str, float]:
    y = np.asarray(y, dtype=int)
    score = 1 / (1 + np.exp(-np.asarray(logit, dtype=float)))
    pred = (score >= 0.5).astype(int)
    tp = float(((pred == 1) & (y == 1)).sum())
    tn = float(((pred == 0) & (y == 0)).sum())
    fp = float(((pred == 1) & (y == 0)).sum())
    fn = float(((pred == 0) & (y == 1)).sum())
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    specificity = tn / (tn + fp) if tn + fp else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    out = {
        "accuracy": float((tp + tn) / max(len(y), 1)),
        "balanced_accuracy": float((recall + specificity) / 2),
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "log_loss": float(-(y * np.log(np.clip(score, 1e-6, 1)) + (1 - y) * np.log(np.clip(1 - score, 1e-6, 1))).mean()),
    }
    out["auc"] = _binary_auc(y, score)
    return out


def multiclass_metrics(y: np.ndarray, logits: np.ndarray, classes: list[str]) -> dict[str, object]:
    y = np.asarray(y, dtype=int)
    probs = softmax(np.asarray(logits, dtype=float), axis=1)
    pred = probs.argmax(axis=1)
    per_f1 = []
    recalls = []
    for idx in range(len(classes)):
        tp = float(((pred == idx) & (y == idx)).sum())
        fp = float(((pred == idx) & (y != idx)).sum())
        fn = float(((pred != idx) & (y == idx)).sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        recalls.append(recall)
        per_f1.append(2 * precision * recall / (precision + recall) if precision + recall else 0.0)
    out: dict[str, object] = {
        "accuracy": float((pred == y).mean()),
        "balanced_accuracy": float(np.mean(recalls)),
        "macro_f1": float(np.mean(per_f1)),
        "per_class_f1": {cls: float(val) for cls, val in zip(classes, per_f1)},
        "log_loss": float(-np.log(np.clip(probs[np.arange(len(y)), y], 1e-6, 1)).mean()),
    }
    aucs = [_binary_auc((y == idx).astype(int), probs[:, idx]) for idx in range(len(classes))]
    out["macro_auc"] = float(np.nanmean(aucs)) if np.isfinite(aucs).any() else float("nan")
    return out


def _binary_auc(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    pos = score[y == 1]
    neg = score[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    wins = 0.0
    for val in pos:
        wins += float((val > neg).sum()) + 0.5 * float((val == neg).sum())
    return float(wins / (len(pos) * len(neg)))


def welch_z_p(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan")
    stat, pvalue = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit")
    return float(stat), float(pvalue)


def correlation_z_p(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan"), float("nan")
    r, pvalue = stats.pearsonr(x, y)
    z = 0.5 * math.log((1 + np.clip(r, -0.999999, 0.999999)) / (1 - np.clip(r, -0.999999, 0.999999)))
    return float(z), float(pvalue)


def dataframe_to_records_safe(df: pd.DataFrame) -> list[dict[str, object]]:
    return df.replace({np.nan: None}).to_dict(orient="records")
