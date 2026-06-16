from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from .abundance_data import AbundanceTable
from .abundance_mil import BagDataset, ScSurvivalMILNet, build_bags_from_counts
from .abundance_result import AbundanceResult


@dataclass
class ScSurvivalDataset:
    table: AbundanceTable
    bags: BagDataset
    outcome: pd.DataFrame
    mode: str

    @property
    def x(self) -> np.ndarray:
        return np.asarray([bag.mean(axis=0) for bag in self.bags.bags], dtype=float)


ScSurvivalModel = ScSurvivalMILNet
ScSurvivalTrainer = None
ScSurvivalResult = AbundanceResult


def make_scsurvival_dataset(
    table: AbundanceTable,
    mode: str,
    outcome: pd.DataFrame,
    covariates: list[str] | None = None,
) -> ScSurvivalDataset:
    return ScSurvivalDataset(
        table=table,
        bags=build_bags_from_counts(table),
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
