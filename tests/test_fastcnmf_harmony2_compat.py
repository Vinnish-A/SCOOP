from __future__ import annotations

import numpy as np

from fastcnmf.harmony2_compat import (
    build_fixed_lamb,
    moe_correct_ridge_batched,
    moe_correct_ridge_fast,
)


def test_fixed_lamb_reconstructs_legacy_penalty_values() -> None:
    lamb = build_fixed_lamb(1, [2, 1])
    assert lamb.shape == (4,)
    np.testing.assert_array_equal(lamb, np.array([0.0, 1.0, 1.0, 1.0]))


def test_batched_moe_matches_cluster_loop_with_diagonal_lambda() -> None:
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


def test_vector_lambda_is_diagonal_and_empty_cluster_contributes_nothing() -> None:
    x = np.array([[2., 1.], [1., 3.], [4., 2.], [1., 4.]])
    phi = np.array([[1., 1., 1., 1.], [1., 1., 0., 0.], [0., 0., 1., 1.]])
    r = np.array([[.7, .3], [.2, .8], [.4, .6], [.8, .2]])
    penalty = np.array([0., 1., 1.])
    # Independent original cNMF MOE equation with Harmony's diagonal lambda.
    signed = x.copy()
    for weights in r.T:
        phi_r = phi * weights
        coeff = np.linalg.solve(phi_r @ phi.T + np.diag(penalty), phi_r @ x)
        coeff[0] = 0
        signed -= phi_r.T @ coeff
    expected = np.maximum(signed, 0)
    np.testing.assert_allclose(moe_correct_ridge_fast(x, r, phi, penalty), expected)
    np.testing.assert_allclose(moe_correct_ridge_batched(x, r, phi, penalty), expected)
    with_empty = np.column_stack([r, np.zeros(len(x))])
    np.testing.assert_allclose(moe_correct_ridge_fast(x, with_empty, phi, penalty), expected)
    np.testing.assert_allclose(moe_correct_ridge_batched(x, with_empty, phi, penalty), expected)
    near_empty = np.column_stack([r, np.full(len(x), 1e-12)])
    np.testing.assert_allclose(moe_correct_ridge_batched(x, near_empty, phi, penalty),
                               moe_correct_ridge_fast(x, near_empty, phi, penalty), rtol=1e-10, atol=1e-10)


def test_all_empty_clusters_match_normal_nonnegative_clipping() -> None:
    x = np.array([[2., -1.], [0., 3.]])
    phi = np.array([[1., 1.], [1., 0.], [0., 1.]])
    empty = np.zeros((2, 2))
    penalty = np.array([0., 1., 1.])
    expected = np.maximum(x, 0.)
    np.testing.assert_array_equal(moe_correct_ridge_fast(x, empty, phi, penalty), expected)
    np.testing.assert_array_equal(moe_correct_ridge_batched(x, empty, phi, penalty), expected)
