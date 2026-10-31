import json

import numpy as np
import pandas as pd
import pytest
import yaml
from anndata import AnnData, read_h5ad
from scipy import sparse

from scoop.atlas import integrate


def _run(tmp_path, scvi=None, split=None):
    rng = np.random.default_rng(0)
    counts = sparse.csr_matrix(rng.poisson(2.0, (60, 20)).astype(np.float32))
    obs = pd.DataFrame({"donor_id": np.repeat(["d1", "d2", "d3"], 20), "batch": np.tile(["b1", "b2"], 30)},
                       index=[f"c{i}" for i in range(60)])
    source = tmp_path / "in.h5ad"
    AnnData(counts, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(20)]), layers={"counts": counts}).write_h5ad(source)
    root = tmp_path / "runs/r1"
    (root / "config").mkdir(parents=True, exist_ok=True)
    config = root / "config/run.yaml"
    config.write_text(yaml.safe_dump({"atlas": {
        "split": {"unit": "donor_id", **(split or {})},
        "scvi": {"n_latent": 2, "max_epochs": 1, "batch_size": 32, "accelerator": "cpu", "devices": 1, **(scvi or {})},
        "graph": {"n_neighbors": 5, "leiden_resolutions": [1.0], "leiden_seeds": [0],
                  "leiden_search": {"backend": "scanpy"}}}}))
    return integrate(config_path=config, input_path=source), root


def test_scvi_batch_key_must_be_declared(tmp_path):
    with pytest.raises(ValueError, match="atlas.scvi.batch_key"):
        _run(tmp_path, scvi={"batch_key": None})


def test_held_out_units_are_projected_not_trained(tmp_path):
    out, root = _run(tmp_path, scvi={"batch_key": "batch"}, split={"held_out_query_units": ["d3"]})
    adata = read_h5ad(out)
    assert adata.obs.is_projected.tolist() == [False] * 40 + [True] * 20
    assert adata.obsm["X_scVI"].shape == (60, 2) and np.isfinite(adata.obsm["X_scVI"]).all()
    params = json.loads((root / "atlas/scvi-run.json").read_text())["parameters"]
    assert params["batch_key"] == "batch" and params["held_out_query_units"] == ["d3"]
    assert (params["n_training_cells"], params["n_projected_cells"]) == (40, 20)


def test_held_out_units_must_exist_and_share_training_batches(tmp_path):
    with pytest.raises(ValueError, match="not found"):
        _run(tmp_path, scvi={"batch_key": "batch"}, split={"held_out_query_units": ["d9"]})
    with pytest.raises(ValueError, match="absent from the training cells"):
        _run(tmp_path, scvi={"batch_key": "donor_id"}, split={"held_out_query_units": ["d3"]})
