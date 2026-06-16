from __future__ import annotations

import numpy as np

from fastde.abundance_loss import cox_partial_nll


def test_cox_loss_finite() -> None:
    loss = cox_partial_nll(np.array([0.2, 0.1, -0.3]), np.array([3.0, 2.0, 1.0]), np.array([1, 0, 1]))
    assert np.isfinite(loss)


def test_cox_loss_lower_for_correct_risk_order() -> None:
    time = np.array([1.0, 2.0, 3.0])
    event = np.array([1, 1, 1])
    correct = cox_partial_nll(np.array([3.0, 2.0, 1.0]), time, event)
    wrong = cox_partial_nll(np.array([1.0, 2.0, 3.0]), time, event)
    assert correct < wrong
