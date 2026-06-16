from __future__ import annotations

import numpy as np
from scipy.special import expit, logsumexp


def cox_partial_nll(risk: np.ndarray, time: np.ndarray, event: np.ndarray) -> float:
    risk = np.asarray(risk, dtype=float).reshape(-1)
    time = np.asarray(time, dtype=float).reshape(-1)
    event = np.asarray(event, dtype=float).reshape(-1)
    if risk.shape[0] != time.shape[0] or risk.shape[0] != event.shape[0]:
        raise ValueError("risk, time, and event must have the same length")
    if event.sum() < 1:
        raise ValueError("Cox partial likelihood requires at least one event")
    order = np.argsort(-time)
    sorted_risk = risk[order]
    sorted_event = event[order]
    log_risk = np.logaddexp.accumulate(sorted_risk)
    observed = sorted_event > 0
    loss = -(sorted_risk[observed] - log_risk[observed]).mean()
    if not np.isfinite(loss):
        raise FloatingPointError("non-finite Cox partial likelihood")
    return float(loss)


def cox_gradient(
    x: np.ndarray,
    beta: np.ndarray,
    time: np.ndarray,
    event: np.ndarray,
    l2: float = 0.0,
) -> tuple[float, np.ndarray]:
    x = np.asarray(x, dtype=float)
    beta = np.asarray(beta, dtype=float)
    order = np.argsort(-np.asarray(time, dtype=float))
    xs = x[order]
    es = np.asarray(event, dtype=float)[order]
    if es.sum() < 1:
        raise ValueError("Cox partial likelihood requires at least one event")
    risk = xs @ beta
    max_prefix = np.maximum.accumulate(risk)
    shifted = np.exp(risk - max_prefix)
    cum_shifted = np.cumsum(shifted)
    # Recompute stable weighted risk-set means with a direct mask for small sample-level data.
    grad = np.zeros_like(beta)
    losses = []
    event_count = int(es.sum())
    for i, is_event in enumerate(es > 0):
        if not is_event:
            continue
        rr = risk[: i + 1]
        weights = np.exp(rr - rr.max())
        denom = weights.sum()
        mean_x = (weights[:, None] * xs[: i + 1]).sum(axis=0) / denom
        grad -= xs[i] - mean_x
        losses.append(-(risk[i] - (np.log(denom) + rr.max())))
    grad = grad / max(event_count, 1) + float(l2) * beta
    loss = float(np.mean(losses) + 0.5 * float(l2) * np.dot(beta, beta))
    return loss, grad


def binary_bce_with_logits(logit: np.ndarray, y: np.ndarray, weight: np.ndarray | None = None) -> float:
    logit = np.asarray(logit, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    loss = np.maximum(logit, 0) - logit * y + np.log1p(np.exp(-np.abs(logit)))
    if weight is not None:
        loss = loss * np.asarray(weight, dtype=float).reshape(-1)
    return float(np.mean(loss))


def multiclass_cross_entropy(logits: np.ndarray, y: np.ndarray) -> float:
    logits = np.asarray(logits, dtype=float)
    y = np.asarray(y, dtype=int).reshape(-1)
    log_probs = logits - logsumexp(logits, axis=1, keepdims=True)
    return float(-log_probs[np.arange(len(y)), y].mean())


def multinomial_nll(counts: np.ndarray, logits: np.ndarray) -> float:
    counts = np.asarray(counts, dtype=float)
    logits = np.asarray(logits, dtype=float)
    log_probs = logits - logsumexp(logits, axis=1, keepdims=True)
    totals = np.maximum(counts.sum(axis=1), 1.0)
    return float((-(counts * log_probs).sum(axis=1) / totals).mean())


def sigmoid(x: np.ndarray) -> np.ndarray:
    return expit(x)
