from __future__ import annotations

import numpy as np

from fastcnmf.harmony2_compat import (
    build_fixed_lamb,
    moe_correct_ridge_batched,
    moe_correct_ridge_fast,
)


def test_fixed_lamb_preserves_cnmf_vector_broadcast_contract() -> None:
    lamb = build_fixed_lamb(1, [2, 1])
    assert lamb.shape == (4,)
    np.testing.assert_array_equal(lamb, np.array([0.0, 1.0, 1.0, 1.0]))


def test_batched_moe_matches_cluster_loop_with_vector_lambda() -> None:
    rng = np.random.default_rng(11)
    x = rng.gamma(shape=1.5, scale=1.0, size=(8, 5))
    r = rng.random((8, 3))
    r /= r.sum(axis=1, keepdims=True)
    batches = np.array([0, 1, 0, 2, 1, 2, 0, 1])
    phi = np.zeros((4, 8))
    phi[0] = 1
    phi[batches + 1, np.arange(8)] = 1
    lamb = np.array([0.0, 1.0, 1.0, 1.0])

    loop = moe_correct_ridge_fast(x, r, phi, lamb)
    batched = moe_correct_ridge_batched(x, r, phi, lamb)

    np.testing.assert_allclose(loop, batched, rtol=1e-10, atol=1e-10)
