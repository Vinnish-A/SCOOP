import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scsp_agent_sop.measurement import FeatureMeasurement, feature_measurement_path
from scsp_agent_sop.storage import register_file

MEASURED = np.array([[1, 1, 0, 1, 1, 1, 1, 1, 0],
                     [1, 0, 0, 1, 1, 1, 1, 1, 1],
                     [1, 1, 1, 1, 1, 1, 1, 1, 1]], bool)


def _adata(samples=("s1", "s1", "s2", "s3")):
    return AnnData(np.zeros((len(samples), MEASURED.shape[1])),
                   obs=pd.DataFrame({"sample": list(samples)}, index=[f"c{i}" for i in range(len(samples))]),
                   var=pd.DataFrame(index=[f"g{i}" for i in range(MEASURED.shape[1])]))


def _measurement(adata):
    reasons = np.where(MEASURED, 1, 4).astype(np.uint8)
    return FeatureMeasurement.from_var("sample", ["s1", "s2", "s3"], adata.var_names, MEASURED,
                                       reasons=reasons, reason_legend={1: "measured", 4: "unmeasured_sample_structural"})


def test_round_trip_keeps_mask_reasons_and_packs_bits(tmp_path):
    adata = _adata()
    path = _measurement(adata).save(tmp_path / "measurement.npz")
    with np.load(path) as z:
        assert z["measured_packbits"].dtype == np.uint8 and z["measured_packbits"].shape == (3, 2)
    loaded = FeatureMeasurement.load(path)
    np.testing.assert_array_equal(loaded.measured, MEASURED)
    assert list(loaded.samples) == ["s1", "s2", "s3"] and loaded.sample_key == "sample"
    assert loaded.reason_legend == {1: "measured", 4: "unmeasured_sample_structural"}
    np.testing.assert_array_equal(loaded.reasons, np.where(MEASURED, 1, 4))
    loaded.check(adata)


def test_var_axis_and_sample_coverage_are_checked():
    adata = _adata()
    fm = _measurement(adata)
    renamed = adata.copy()
    renamed.var_names = [f"x{i}" for i in range(adata.n_vars)]
    with pytest.raises(ValueError, match="var_sha256"):
        fm.check(renamed)
    with pytest.raises(ValueError, match="does not cover 1 cells"):
        fm.check(_adata(("s1", "s4")))


def test_measured_views_follow_sample_codes():
    adata = _adata()
    fm = _measurement(adata)
    samples = adata.obs["sample"]
    np.testing.assert_array_equal(fm.measured_in_all(samples), MEASURED.all(0))
    np.testing.assert_array_equal(fm.measured_in_all(["s3", "s3"]), np.ones(9, bool))
    np.testing.assert_array_equal(fm.measured_cells(samples), [4, 3, 1, 4, 4, 4, 4, 4, 2])
    np.testing.assert_array_equal(fm.cell_mask(samples, [1, 8]), [[1, 0], [1, 0], [0, 1], [1, 1]])


def test_discovery_prefers_explicit_path_then_registry(tmp_path):
    adata = _adata()
    assert feature_measurement_path(adata) is None
    path = _measurement(adata).save(tmp_path / "fm.npz")
    register_file(adata, key="feature_measurement", path=path, category="artifacts")
    assert feature_measurement_path(adata) == path
    assert feature_measurement_path(adata, tmp_path / "other.npz") == tmp_path / "other.npz"
