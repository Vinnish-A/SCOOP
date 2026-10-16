from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse
from sklearn.metrics import adjusted_rand_score

from scsp_agent_sop.parallel_leiden import run_jobs


def test_weighted_components_parallel_and_resume(tmp_path):
    pytest.importorskip("networkit")
    blocks = [np.ones((20, 20)) - np.eye(20) for _ in range(3)]
    graph = sparse.block_diag(blocks, format="csr")
    cfg = {"workers": 2, "threads_per_worker": 2, "iterations": 5, "checkpoint_dir": str(tmp_path)}
    jobs = [(0.5, 0), (1.0, 1)]
    results = {(r, s): labels for r, s, labels in run_jobs(graph, jobs, cfg)}
    for labels in results.values():
        assert adjusted_rand_score(np.repeat(np.arange(3), 20), labels) == 1.0
    resumed = {(r, s): labels for r, s, labels in run_jobs(graph, jobs, cfg)}
    for key in results:
        np.testing.assert_array_equal(results[key], resumed[key])
    with pytest.raises(ValueError, match="Checkpoint settings changed"):
        list(run_jobs(graph, jobs, dict(cfg, iterations=6)))


def test_directed_graph_rejected(tmp_path):
    pytest.importorskip("networkit")
    graph = sparse.csr_matrix([[0, 1], [0, 0]])
    with pytest.raises(ValueError, match="symmetric"):
        list(run_jobs(graph, [(1., 0)], {"checkpoint_dir": str(tmp_path)}))
