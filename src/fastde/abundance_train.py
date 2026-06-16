from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .abundance_model import ScSurvivalLikeModel


@dataclass
class AbundanceTrainingConfig:
    learning_rate: float = 1e-3
    weight_decay: float = 1e-4
    max_epochs: int = 500
    random_seed: int = 0
    binary_class_weight: str | None = "balanced"


@dataclass
class AbundanceTrainer:
    mode: str
    config: AbundanceTrainingConfig

    def fit(
        self,
        x: pd.DataFrame,
        y: np.ndarray | None = None,
        time: np.ndarray | None = None,
        event: np.ndarray | None = None,
    ) -> ScSurvivalLikeModel:
        model = ScSurvivalLikeModel(
            mode=self.mode,
            learning_rate=self.config.learning_rate,
            weight_decay=self.config.weight_decay,
            max_epochs=self.config.max_epochs,
            random_seed=self.config.random_seed,
            binary_class_weight=self.config.binary_class_weight,
        )
        values = x.to_numpy(dtype=float)
        if self.mode in {"binary", "condition"}:
            if y is None:
                raise ValueError("binary training requires labels")
            return model.fit_binary(values, y)
        if self.mode == "multiclass":
            if y is None:
                raise ValueError("multiclass training requires labels")
            return model.fit_multiclass(values, y)
        if self.mode == "continuous":
            if y is None:
                raise ValueError("continuous training requires labels")
            return model.fit_continuous(values, y)
        if self.mode == "survival":
            if time is None or event is None:
                raise ValueError("survival training requires time and event")
            return model.fit_survival(values, time, event)
        raise ValueError(f"unknown abundance training mode: {self.mode}")
