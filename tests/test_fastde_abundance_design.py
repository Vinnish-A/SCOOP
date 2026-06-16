from __future__ import annotations

import pandas as pd

from fastde.abundance_data import AbundanceTable
from fastde.abundance_design import build_feature_design, encode_multiclass_labels


def test_multiclass_design_reference_level() -> None:
    counts = pd.DataFrame([[5, 1], [1, 5], [4, 2]], index=["S1", "S2", "S3"], columns=["A", "B"])
    meta = pd.DataFrame({"subtype": ["case", "control", "case"], "batch": ["x", "x", "y"]}, index=counts.index)
    table = AbundanceTable(counts, counts.div(counts.sum(axis=1), axis=0), meta)
    y, classes, keep, ref = encode_multiclass_labels(meta, "subtype", "control")
    design = build_feature_design(table, covariates=["batch"])
    assert classes[ref] == "control"
    assert list(keep) == ["S1", "S2", "S3"]
    assert design.features.shape[0] == 3
    assert "batch_y" in design.features.columns
