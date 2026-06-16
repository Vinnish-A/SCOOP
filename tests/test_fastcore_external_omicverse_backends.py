from __future__ import annotations

import sys
import types

import numpy as np
from anndata import AnnData
from scipy import sparse

from fastcore.backend_plan import FastCorePlan
from fastcore.backends.omicverse_mixed import run_omicverse_mixed_core
from scsp_agent_sop.core_runner import run_core_pipeline
from scsp_agent_sop.storage import init_file_registry


def _adata() -> AnnData:
    adata = AnnData(sparse.csr_matrix(np.arange(60, dtype=np.float32).reshape(10, 6) + 1))
    adata.obs_names = [f"cell{i}" for i in range(adata.n_obs)]
    adata.var_names = [f"gene{i}" for i in range(adata.n_vars)]
    adata.layers["counts"] = adata.X.copy()
    init_file_registry(adata, "test")
    return adata


def _cfg() -> dict:
    return {
        "run": {"run_id": "test", "random_seed": 0},
        "keys": {"batch_candidates": []},
        "qc": {"counts_layer": "counts"},
        "core": {
            "engine": "fastcore",
            "fallback_engine": "scanpy_legacy",
            "batch_correction": {"method": "harmony2", "max_iter_harmony": 20},
            "fastcore": {
                "fallback_backend": "scanpy_legacy",
                "omicverse": {
                    "preprocess_mode": "shiftlog",
                    "target_sum": 500000,
                    "n_hvgs": 4,
                    "n_pcs": 3,
                    "neighbors": {"n_neighbors": 3, "n_pcs": 3, "method": "auto", "mixed_transformer": "pyg"},
                    "umap": {"min_dist": 0.3, "method": "auto"},
                    "leiden": {"resolutions": [0.6], "default_resolution": 0.6},
                    "rust_oom": {"preprocess_mode": "shiftlog|pearson"},
                },
            },
        },
    }


class _FakeSettings:
    def __init__(self, calls: list[str]):
        self.calls = calls
        self.mode = "cpu"

    def cpu_init(self):
        self.mode = "cpu"
        self.calls.append("cpu_init")

    def cpu_gpu_mixed_init(self):
        self.mode = "cpu-gpu-mixed"
        self.calls.append("mixed_init")

class _FakePP:
    def __init__(self, calls: list[str]):
        self.calls = calls

    def preprocess(self, adata, **kwargs):
        self.calls.append(f"preprocess:{kwargs['mode']}")
        adata.var["highly_variable"] = [True, True, True, True, False, False]
        adata.var["means"] = np.arange(adata.n_vars, dtype=float)
        adata.layers["lognorm"] = adata.X.copy()

    def scale(self, adata, layers_add="scaled", **kwargs):
        self.calls.append(f"scale:{layers_add}:{kwargs.get('use_implicit_centering')}")
        adata.layers[layers_add] = adata.X.copy()
        if kwargs.get("use_implicit_centering"):
            adata.uns["_scaled_implicit"] = object()

    def pca(self, adata, n_pcs=3, **kwargs):
        self.calls.append(f"pca:{n_pcs}")
        adata.obsm["X_pca"] = np.arange(adata.n_obs * n_pcs, dtype=np.float32).reshape(adata.n_obs, n_pcs)
        adata.uns["pca"] = {"variance_ratio": np.ones(n_pcs), "variance": np.ones(n_pcs)}

    def neighbors(self, adata, key_added="neighbors_identity", **kwargs):
        self.calls.append(f"neighbors:{kwargs.get('method')}:{kwargs.get('transformer')}")
        graph = sparse.eye(adata.n_obs, format="csr")
        adata.obsp[f"{key_added}_connectivities"] = graph
        adata.obsp[f"{key_added}_distances"] = graph

    def umap(self, adata, **kwargs):
        self.calls.append(f"umap:{kwargs.get('method')}")
        adata.obsm["X_umap"] = np.zeros((adata.n_obs, 2), dtype=np.float32)

    def leiden(self, adata, key_added="cluster_identity", **kwargs):
        self.calls.append(f"leiden:{kwargs['resolution']}")
        adata.obs[key_added] = ["0"] * (adata.n_obs // 2) + ["1"] * (adata.n_obs - adata.n_obs // 2)


class _FakeOOM:
    def __init__(self, adata: AnnData, calls: list[str]):
        self._adata = adata
        self._calls = calls

    def __getattr__(self, name):
        return getattr(self._adata, name)

    def to_adata(self):
        self._calls.append("to_adata")
        return self._adata

    def close(self):
        self._calls.append("close")


def _install_fake_omicverse(monkeypatch, *, rust_input: AnnData | None = None) -> list[str]:
    calls: list[str] = []
    fake = types.SimpleNamespace()
    fake.settings = _FakeSettings(calls)
    fake.pp = _FakePP(calls)
    fake.set_seed = lambda seed, verbose=False: calls.append(f"seed:{seed}")

    def read(path, backend="python", **kwargs):
        calls.append(f"read:{backend}:{path}")
        return _FakeOOM(rust_input or _adata(), calls)

    fake.read = read
    monkeypatch.setitem(sys.modules, "omicverse", fake)
    return calls


def test_external_mixed_backend_runs_omicverse_steps(monkeypatch, tmp_path):
    calls = _install_fake_omicverse(monkeypatch)
    adata = _adata()

    result = run_omicverse_mixed_core(adata, _cfg(), tmp_path)

    assert result["backend"] == "fastcore_mixed"
    assert calls[:2] == ["seed:0", "mixed_init"]
    assert "X_pca_harmony_identity" in adata.obsm
    assert "connectivities_identity" in adata.obsp
    assert "cluster_identity" in adata.obs
    assert "_scaled_implicit" not in adata.uns
    assert any(call == "neighbors:None:pyg" for call in calls)


def test_rust_oom_backend_is_path_based_and_writes_output(monkeypatch, tmp_path):
    calls = _install_fake_omicverse(monkeypatch, rust_input=_adata())

    def fake_plan(cfg, *, adata=None, input_path=None, capabilities=None):
        return FastCorePlan("fastcore_oom", False, "scanpy_legacy", [], {})

    monkeypatch.setattr("scsp_agent_sop.core_runner.plan_fastcore_backend", fake_plan)
    output = tmp_path / "adata_core.h5ad"

    result = run_core_pipeline(None, _cfg(), tmp_path, input_path=tmp_path / "input.h5ad", output_path=output)

    assert result["backend"] == "fastcore_oom"
    assert output.exists()
    assert calls[0].startswith("read:rust:")
    assert "to_adata" in calls
    assert calls[-1] == "close"
