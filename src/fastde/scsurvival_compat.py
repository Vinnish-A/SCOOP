from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from .abundance_data import AbundanceTable
from .abundance_design import FeatureDesign, build_feature_design
from .abundance_model import ScSurvivalLikeModel
from .abundance_result import AbundanceResult
from .abundance_train import AbundanceTrainer


@dataclass
class ScSurvivalDataset:
    table: AbundanceTable
    design: FeatureDesign
    outcome: pd.DataFrame
    mode: str

    @property
    def x(self) -> np.ndarray:
        return self.design.features.to_numpy(dtype=float)


ScSurvivalModel = ScSurvivalLikeModel
ScSurvivalTrainer = AbundanceTrainer
ScSurvivalResult = AbundanceResult


def make_scsurvival_dataset(
    table: AbundanceTable,
    mode: str,
    outcome: pd.DataFrame,
    transform: str = "clr",
    pseudocount: float = 0.5,
    covariates: list[str] | None = None,
) -> ScSurvivalDataset:
    return ScSurvivalDataset(
        table=table,
        design=build_feature_design(table, transform=transform, pseudocount=pseudocount, covariates=covariates),
        outcome=outcome,
        mode=mode,
    )


def to_fastde_result(result: ScSurvivalResult) -> dict[str, Any]:
    return {
        "mode": result.mode,
        "results": result.results,
        "predictions": result.predictions,
        "metrics": result.metrics,
        "manifest": result.manifest,
    }
