from __future__ import annotations

import numpy as np
from scipy import sparse

from fastcore.quality import graph_density, pc_subspace_cosine
from scoop_fast.registry import get_engine, list_engines


def test_fastcore_registry_exposes_expected_engines():
    engine_ids = {engine.engine_id for engine in list_engines()}
    assert {
        "fastcore.preprocess",
        "fastcore.core_pipeline",
        "fastcore.graph_embed_cluster",
        "fastcore.quality_compare",
    }.issubset(engine_ids)
    assert get_engine("fastcore.core_pipeline").task_type == "core_pipeline"


def test_pc_subspace_cosine_is_sign_invariant():
    x = np.eye(4, 2)
    y = x.copy()
    y[:, 0] *= -1
    assert pc_subspace_cosine(x, y) == 1.0


def test_graph_density_sparse():
    graph = sparse.csr_matrix(
        np.array(
            [
                [0, 1, 0],
                [1, 0, 1],
                [0, 1, 0],
            ]
        )
    )
    assert graph_density(graph) == 4 / 6
