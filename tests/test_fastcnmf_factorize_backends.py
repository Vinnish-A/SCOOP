from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml

from fastcnmf.fast_factorize import run_fast_factorize


def _save_df_to_npz(obj: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, data=obj.values, index=obj.index.values, columns=obj.columns.values)


def test_minibatch_factorize_writes_cnmf_spectra(tmp_path: Path) -> None:
    pytest.importorskip("sklearn")

    output_dir = tmp_path / "nmf"
    run_name = "tiny"
    tmp_dir = output_dir / run_name / "cnmf_tmp"
    tmp_dir.mkdir(parents=True)

    x = np.abs(np.random.default_rng(1).normal(size=(20, 8))).astype(np.float64)
    ad.AnnData(
        X=x,
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(x.shape[0])]),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(x.shape[1])]),
    ).write_h5ad(tmp_dir / f"{run_name}.norm_counts.h5ad")
    params = pd.DataFrame(
        [[3, 0, 11, False]],
        columns=["n_components", "iter", "nmf_seed", "completed"],
    )
    _save_df_to_npz(params, tmp_dir / f"{run_name}.nmf_params.df.npz")
    with (tmp_dir / f"{run_name}.nmf_idvrun_params.yaml").open("w", encoding="utf-8") as handle:
        yaml.dump(
            {
                "alpha_W": 0.0,
                "alpha_H": 0.0,
                "l1_ratio": 0.0,
                "beta_loss": "frobenius",
                "solver": "cd",
                "tol": 1e-4,
                "max_iter": 10,
                "init": "random",
            },
            handle,
        )

    result = run_fast_factorize(
        output_dir=output_dir,
        run_name=run_name,
        workers=1,
        backend="minibatch",
        minibatch_batch_size=8,
    )

    output = tmp_dir / f"{run_name}.spectra.k_3.iter_0.df.npz"
    assert result["backend"] == "minibatch"
    assert output.exists()
    with np.load(output, allow_pickle=True) as payload:
        spectra = pd.DataFrame(**payload)
    assert spectra.shape == (3, 8)
    assert list(spectra.columns) == [f"gene{i}" for i in range(8)]
