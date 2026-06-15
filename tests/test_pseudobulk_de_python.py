from __future__ import annotations

import numpy as np
import pandas as pd

from scsp_agent_sop.de import benjamini_hochberg, compare_de_tables, run_pseudobulk_logcpm_welch


def test_benjamini_hochberg_monotone_adjustment() -> None:
    p = np.array([0.03, 0.001, 0.02, 0.5])
    fdr = benjamini_hochberg(p)
    assert np.all((fdr >= 0) & (fdr <= 1))
    ordered = fdr[np.argsort(p)]
    assert np.all(np.diff(ordered) >= -1e-12)


def test_run_pseudobulk_logcpm_welch_detects_direction() -> None:
    counts = pd.DataFrame(
        {
            "UP": [100, 110, 105, 420, 440, 430],
            "FLAT": [100, 100, 105, 100, 98, 102],
            "LOW": [20, 18, 21, 19, 20, 18],
        },
        index=[f"s{i}" for i in range(6)],
    )
    meta = pd.DataFrame({"condition": ["ctrl", "ctrl", "ctrl", "test", "test", "test"]}, index=counts.index)

    result = run_pseudobulk_logcpm_welch(counts, meta, condition_col="condition", ctrl_group="ctrl", test_group="test", min_count=1)
    table = result.table.set_index("gene")

    assert {"gene", "logFC", "logCPM", "PValue", "FDR"}.issubset(result.table.columns)
    assert table.loc["UP", "logFC"] > 0
    assert np.isfinite(table.loc["UP", "PValue"])
    assert np.isfinite(table.loc["UP", "FDR"])
    assert result.manifest["method"] == "python_logcpm_welch"


def test_compare_de_tables_reports_overlap() -> None:
    left = pd.DataFrame({"gene": ["A", "B", "C"], "logFC": [1.0, -1.0, 0.2], "PValue": [0.001, 0.02, 0.5]})
    right = pd.DataFrame({"gene": ["A", "B", "C"], "logFC": [0.8, -0.5, 0.1], "PValue": [0.002, 0.03, 0.4]})

    report = compare_de_tables(left, right, top_n=2)

    assert report["n_common_genes"] == 3
    assert report["top_n_overlap_fraction"] == 1.0
    assert report["logfc_sign_agreement"] == 1.0
