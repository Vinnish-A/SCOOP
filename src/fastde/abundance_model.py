from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .abundance_loss import cox_gradient


@dataclass
class ScSurvivalLikeModel:
    mode: str
    learning_rate: float = 1e-3
    weight_decay: float = 1e-4
    max_epochs: int = 500
    random_seed: int = 0
    binary_class_weight: str | dict[int, float] | None = "balanced"
    estimator: object | None = None
    coef_: np.ndarray | None = None
    intercept_: np.ndarray | float | None = None
    history_: list[dict[str, float]] = field(default_factory=list)
    mean_: np.ndarray | None = None
    scale_: np.ndarray | None = None

    def _scale_fit(self, x: np.ndarray) -> np.ndarray:
        arr = np.asarray(x, dtype=float)
        self.mean_ = arr.mean(axis=0)
        self.scale_ = arr.std(axis=0)
        self.scale_[self.scale_ == 0] = 1.0
        return (arr - self.mean_) / self.scale_

    def _scale(self, x: np.ndarray) -> np.ndarray:
        if self.mean_ is None or self.scale_ is None:
            raise RuntimeError("model has not been fitted")
        return (np.asarray(x, dtype=float) - self.mean_) / self.scale_

    def fit_binary(self, x: np.ndarray, y: np.ndarray) -> "ScSurvivalLikeModel":
        xs = self._scale_fit(x)
        y = np.asarray(y, dtype=float)
        rng = np.random.default_rng(self.random_seed)
        beta = rng.normal(0, 1e-3, size=xs.shape[1])
        intercept = 0.0
        if self.binary_class_weight == "balanced":
            pos_w = len(y) / max(2 * y.sum(), 1)
            neg_w = len(y) / max(2 * (len(y) - y.sum()), 1)
            weights = np.where(y == 1, pos_w, neg_w)
        else:
            weights = np.ones_like(y)
        for epoch in range(int(self.max_epochs)):
            logit = xs @ beta + intercept
            prob = 1 / (1 + np.exp(-np.clip(logit, -40, 40)))
            err = (prob - y) * weights
            grad = xs.T @ err / len(y) + self.weight_decay * beta
            grad_i = float(err.mean())
            beta -= self.learning_rate * grad
            intercept -= self.learning_rate * grad_i
            if epoch % 25 == 0:
                loss = np.mean(weights * (np.maximum(logit, 0) - logit * y + np.log1p(np.exp(-np.abs(logit)))))
                self.history_.append({"epoch": float(epoch), "loss": float(loss)})
        self.coef_ = beta
        self.intercept_ = intercept
        return self

    def fit_multiclass(self, x: np.ndarray, y: np.ndarray) -> "ScSurvivalLikeModel":
        xs = self._scale_fit(x)
        y = np.asarray(y, dtype=int)
        classes = np.unique(y)
        n_classes = int(classes.max()) + 1
        target = np.zeros((len(y), n_classes), dtype=float)
        target[np.arange(len(y)), y] = 1.0
        rng = np.random.default_rng(self.random_seed)
        beta = rng.normal(0, 1e-3, size=(n_classes, xs.shape[1]))
        intercept = np.zeros(n_classes, dtype=float)
        for epoch in range(int(self.max_epochs)):
            logits = xs @ beta.T + intercept
            logits = logits - logits.max(axis=1, keepdims=True)
            exp_logits = np.exp(logits)
            probs = exp_logits / exp_logits.sum(axis=1, keepdims=True)
            err = (probs - target) / len(y)
            grad = err.T @ xs + self.weight_decay * beta
            grad_i = err.sum(axis=0)
            beta -= self.learning_rate * grad
            intercept -= self.learning_rate * grad_i
            if epoch % 25 == 0:
                loss = -np.log(np.clip(probs[np.arange(len(y)), y], 1e-8, 1)).mean()
                self.history_.append({"epoch": float(epoch), "loss": float(loss)})
        self.coef_ = beta
        self.intercept_ = intercept
        return self

    def fit_continuous(self, x: np.ndarray, y: np.ndarray) -> "ScSurvivalLikeModel":
        xs = self._scale_fit(x)
        y = np.asarray(y, dtype=float)
        y_mean = float(y.mean())
        yc = y - y_mean
        xtx = xs.T @ xs + max(self.weight_decay, 1e-8) * np.eye(xs.shape[1])
        self.coef_ = np.linalg.solve(xtx, xs.T @ yc)
        self.intercept_ = y_mean
        return self

    def fit_survival(self, x: np.ndarray, time: np.ndarray, event: np.ndarray) -> "ScSurvivalLikeModel":
        xs = self._scale_fit(x)
        rng = np.random.default_rng(self.random_seed)
        beta = rng.normal(0.0, 1e-3, size=xs.shape[1])
        lr = float(self.learning_rate)
        best = float("inf")
        bad_epochs = 0
        patience = max(25, min(100, self.max_epochs // 5))
        for epoch in range(int(self.max_epochs)):
            loss, grad = cox_gradient(xs, beta, time, event, l2=self.weight_decay)
            beta -= lr * grad
            self.history_.append({"epoch": float(epoch), "loss": float(loss)})
            if loss + 1e-8 < best:
                best = loss
                bad_epochs = 0
            else:
                bad_epochs += 1
            if bad_epochs >= patience:
                break
        self.estimator = None
        self.coef_ = beta
        self.intercept_ = 0.0
        return self

    def predict_score(self, x: np.ndarray) -> np.ndarray:
        xs = self._scale(x)
        if self.mode in {"binary", "condition"}:
            assert self.coef_ is not None
            return xs @ self.coef_ + float(self.intercept_ or 0.0)
        if self.mode == "survival":
            assert self.coef_ is not None
            return xs @ self.coef_
        if self.mode == "continuous":
            assert self.coef_ is not None
            return xs @ self.coef_ + float(self.intercept_ or 0.0)
        raise ValueError(f"predict_score is not valid for mode {self.mode}")

    def predict_logits(self, x: np.ndarray) -> np.ndarray:
        xs = self._scale(x)
        if self.mode == "multiclass":
            assert self.coef_ is not None
            return xs @ self.coef_.T + np.asarray(self.intercept_)
        return self.predict_score(x).reshape(-1, 1)
