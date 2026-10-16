from types import SimpleNamespace

import numpy as np
import pytest

from scsp_agent_sop.rctd_adapter import selected_weights


def test_multi_uses_sparse_solution_not_initial_full_fit():
    result = SimpleNamespace(
        weights=np.full((2, 3), 0.3),
        sub_weights=np.array([[0.8, 0.2], [1.0, 0.0]]),
        cell_type_indices=np.array([[2, 0], [1, -1]]),
    )
    np.testing.assert_allclose(selected_weights(result, "multi"), [[0.2, 0.0, 0.8], [0.0, 1.0, 0.0]])
    np.testing.assert_array_equal(selected_weights(result, "full"), result.weights)


def test_doublet_preserves_upstream_pair_weights_separately_from_diagnostic_class():
    result = SimpleNamespace(
        weights=np.full((2, 3), 0.3),
        first_type=np.array([0, 1]),
        second_type=np.array([2, 1]),
        weights_doublet=np.array([[0.7, 0.3], [0.6, 0.4]]),
    )
    np.testing.assert_allclose(selected_weights(result, "doublet"), [[0.7, 0.0, 0.3], [0.0, 1.0, 0.0]])
    with pytest.raises(ValueError):
        selected_weights(result, "unknown")
