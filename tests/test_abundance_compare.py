import numpy as np
import pandas as pd
import pytest
import torch

from fastde.abundance_compare import BinaryCountBagMIL, binary_cv, donor_tests, validate_inputs, run_comparison


def test_donor_only_comparison_does_not_fit_mil(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Donor-only analysis must not fit a classifier")
    monkeypatch.setattr("fastde.abundance_compare.binary_cv", forbidden)
    samples = [f"s{i}" for i in range(8)]
    counts = pd.DataFrame([[6, 4]] * 4 + [[2, 8]] * 4, index=samples, columns=["high", "other"])
    meta = pd.DataFrame({"donor_id": samples, "ACPA": ["Positive"] * 4 + ["Negative"] * 4}, index=samples)
    counts.to_csv(tmp_path / "counts.tsv", sep="\t")
    meta.to_csv(tmp_path / "metadata.tsv", sep="\t")
    result = run_comparison({"run_mil": False, "bootstrap": 100, "seed": 0,
        "cohorts": [{"name": "test", "counts": str(tmp_path / "counts.tsv"),
                     "metadata": str(tmp_path / "metadata.tsv"), "target": "high"}]}, tmp_path / "result")
    assert result.delta_mean.iloc[0] == pytest.approx(.4)
    assert result.donor_permutation_p.iloc[0] == pytest.approx(2 / 70)
    assert "donor_permutation_p_bh" in result
    assert "oof_auc" not in result


def test_compressed_bag_equals_explicit_cell_attention():
    torch.manual_seed(2)
    model = BinaryCountBagMIL(3).eval()
    counts = torch.tensor([[3., 0., 8.], [1., 2., 4.]])
    explicit = []
    for row in counts:
        instances = torch.repeat_interleave(model.eye, row.long(), dim=0)
        h = model.encoder(instances)
        weight = torch.softmax(model.attention(h), dim=0)
        explicit.append(model.head((weight * h).sum(dim=0)).flatten())
    torch.testing.assert_close(model(counts), torch.cat(explicit))
    torch.testing.assert_close(model(counts), model(counts * 5))


def test_duplicate_donors_rejected():
    counts = pd.DataFrame([[1, 9], [2, 8]], index=["a", "b"])
    meta = pd.DataFrame({"donor_id": ["d1", "d1"], "ACPA": ["Positive", "Negative"]}, index=counts.index)
    with pytest.raises(ValueError, match="one bag per donor"):
        validate_inputs(counts, meta)


def test_patient_tests_and_cross_validation_separate_donors():
    x = np.array([[80, 20]] * 6 + [[20, 80]] * 6, dtype=float)
    y = np.array([1] * 6 + [0] * 6)
    result = donor_tests(x, y, 0, bootstrap=100, seed=4)
    assert result["delta_mean"] == pytest.approx(.6)
    assert result["roe_positive"] == pytest.approx(1.6)
    assert result["donor_permutation_p"] == pytest.approx(2 / 924)
    p, removed, metrics, folds = binary_cv(x, y, 0, seeds=[0], folds=3, epochs=100)
    assert metrics["oof_auc"] == 1
    assert np.isfinite(removed).all()
    assert np.mean(p[:6]) > np.mean(p[6:])
    assert sorted(i for f in folds for i in f["test"]) == list(range(12))
    for fold in folds:
        assert not set(fold["train"]) & set(fold["test"])
