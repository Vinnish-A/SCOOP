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


def _numpy_cupy(monkeypatch: pytest.MonkeyPatch) -> None:
    """Run the cupy-mu backend's own code on the CPU."""

    from types import SimpleNamespace
    import sys

    fake = SimpleNamespace(
        float32=np.float32,
        float64=np.float64,
        asarray=np.asarray,
        maximum=np.maximum,
        vdot=np.vdot,
        asnumpy=np.asarray,
        cuda=SimpleNamespace(Stream=SimpleNamespace(null=SimpleNamespace(synchronize=lambda: None))),
        get_default_memory_pool=lambda: SimpleNamespace(free_all_blocks=lambda: None),
    )
    monkeypatch.setitem(sys.modules, "cupy", fake)


def _prepare_tiny(tmp_path: Path, *, norm_dtype: str, max_nmf_iter: int = 50) -> tuple[Path, str]:
    from fastcnmf.fast_prepare import run_fast_prepare

    rng = np.random.default_rng(0)
    x = rng.poisson(rng.gamma(1.0, 1.0, (120, 3)) @ rng.gamma(0.5, 2.0, (3, 40)) * 5).astype(np.float64)
    x[x.sum(axis=1) == 0, 0] = 1
    obs = pd.DataFrame(index=[f"c{i}" for i in range(x.shape[0])])
    var = pd.DataFrame(index=[f"g{i}" for i in range(x.shape[1])])
    ad.AnnData(x, obs=obs, var=var).write_h5ad(tmp_path / "corrected.h5ad")
    ad.AnnData(x / x.sum(axis=1, keepdims=True) * 1e4, obs=obs, var=var).write_h5ad(tmp_path / "tpm.h5ad")
    (tmp_path / "hvgs.txt").write_text("\n".join(var.index), encoding="utf-8")
    output_dir = tmp_path / "nmf"
    run_fast_prepare(
        corrected_h5ad=tmp_path / "corrected.h5ad",
        tpm_h5ad=tmp_path / "tpm.h5ad",
        hvgs_txt=tmp_path / "hvgs.txt",
        output_dir=output_dir,
        run_name="tiny",
        k_values=(3,),
        n_iter=4,
        norm_dtype=norm_dtype,
        max_nmf_iter=max_nmf_iter,
    )
    return output_dir, "tiny"


@pytest.mark.parametrize(("norm_dtype", "mu_dtype"), [("float32", "float32"), ("float64", "float32")])
def test_cupy_mu_spectra_keep_norm_counts_dtype_through_consensus(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, norm_dtype: str, mu_dtype: str
) -> None:
    pytest.importorskip("cnmf")
    from fastcnmf.fast_consensus import run_fast_consensus

    _numpy_cupy(monkeypatch)
    output_dir, run_name = _prepare_tiny(tmp_path, norm_dtype=norm_dtype)
    run_fast_factorize(output_dir=output_dir, run_name=run_name, workers=1, backend="cupy-mu", mu_dtype=mu_dtype)

    spectra_path = output_dir / run_name / "cnmf_tmp" / f"{run_name}.spectra.k_3.iter_0.df.npz"
    with np.load(spectra_path, allow_pickle=True) as payload:
        assert payload["data"].dtype == np.dtype(norm_dtype)
    result = run_fast_consensus(
        output_dir=output_dir, run_name=run_name, workers=1, k_values=(3,), k_selection=False, build_reference=False
    )
    assert result["consensus_results"] == [{"k": 3, "stage": "consensus"}]


def test_lite_consensus_refits_float64_spectra_against_float32_counts(tmp_path: Path) -> None:
    pytest.importorskip("cnmf")
    from fastcnmf.fast_consensus import run_fast_consensus

    output_dir, run_name = _prepare_tiny(tmp_path, norm_dtype="float32")
    run_fast_factorize(output_dir=output_dir, run_name=run_name, workers=1, backend="exact")
    # Spectra written by the pre-fix cupy-mu backend were float64 regardless of X.
    for path in (output_dir / run_name / "cnmf_tmp").glob(f"{run_name}.spectra.k_3.iter_*.df.npz"):
        with np.load(path, allow_pickle=True) as payload:
            frame = pd.DataFrame(**payload)
        _save_df_to_npz(frame.astype(np.float64), path)

    result = run_fast_consensus(
        output_dir=output_dir, run_name=run_name, workers=1, k_values=(3,), k_selection=False, lite=True
    )
    assert result["consensus_results"] == [{"k": 3, "stage": "consensus-lite"}]


def test_cupy_mu_stops_on_tolerance_and_reports_actual_iterations(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _numpy_cupy(monkeypatch)
    output_dir, run_name = _prepare_tiny(tmp_path, norm_dtype="float64", max_nmf_iter=2000)

    converged = run_fast_factorize(output_dir=output_dir, run_name=run_name, workers=1, backend="cupy-mu", mu_dtype="float64")
    params_path = output_dir / run_name / "cnmf_tmp" / f"{run_name}.nmf_idvrun_params.yaml"
    params = yaml.safe_load(params_path.read_text(encoding="utf-8"))
    params.update(tol=0.0, max_iter=30)
    params_path.write_text(yaml.safe_dump(params), encoding="utf-8")
    fixed = run_fast_factorize(output_dir=output_dir, run_name=run_name, workers=1, backend="cupy-mu", mu_dtype="float64")

    for item in converged["results"]:
        assert item["n_iter"] % 10 == 0 and 0 < item["n_iter"] < 2000
        assert np.isfinite(item["reconstruction_err"]) and item["reconstruction_err"] > 0
    assert [item["n_iter"] for item in fixed["results"]] == [30] * 4


def test_exact_factorize_reports_reconstruction_error(tmp_path: Path) -> None:
    from sklearn.decomposition import non_negative_factorization

    output_dir, run_name = _prepare_tiny(tmp_path, norm_dtype="float64")
    result = run_fast_factorize(output_dir=output_dir, run_name=run_name, workers=1, backend="exact")

    first = result["results"][0]
    tmp_dir = output_dir / run_name / "cnmf_tmp"
    x = ad.read_h5ad(tmp_dir / f"{run_name}.norm_counts.h5ad").X
    params = yaml.safe_load((tmp_dir / f"{run_name}.nmf_idvrun_params.yaml").read_text(encoding="utf-8"))
    w, h, _ = non_negative_factorization(x, n_components=3, random_state=first["seed"], **params)
    assert np.isclose(first["reconstruction_err"], np.linalg.norm(x - w @ h), rtol=1e-6)
