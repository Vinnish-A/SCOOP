from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from fastcnmf.preprocess_quality import compare_preprocess_outputs


def test_compare_preprocess_outputs_chunked(tmp_path: Path) -> None:
    obs = pd.DataFrame(index=[f"cell{i}" for i in range(5)])
    var = pd.DataFrame(index=[f"gene{i}" for i in range(4)])
    x = np.arange(20, dtype=float).reshape(5, 4)
    ref = tmp_path / "ref.h5ad"
    cand = tmp_path / "cand.h5ad"
    ad.AnnData(x, obs=obs, var=var).write_h5ad(ref)
    ad.AnnData(x + 0.001, obs=obs, var=var).write_h5ad(cand)

    result = compare_preprocess_outputs(
        reference_h5ad=ref,
        candidate_h5ad=cand,
        output_json=tmp_path / "compare.json",
        chunk_size=2,
    )

    assert result["common_obs"] == 5
    assert result["common_vars"] == 4
    assert result["matrix"]["cosine"] > 0.999
    assert result["matrix"]["pearson"] > 0.999
    assert result["passes_95pct_input_gate"]
