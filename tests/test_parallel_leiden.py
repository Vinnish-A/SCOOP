from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse
from sklearn.metrics import adjusted_rand_score

from scsp_agent_sop.parallel_leiden import run_jobs
from scsp_agent_sop import parallel_leiden


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_initialize_preserves_nonuniform_weights_duplicates_and_isolate(tmp_path, dtype):
    pytest.importorskip("networkit")
    # Duplicate (0,1) values sum to0.75. Diagonal99 is ignored by the existing
    # upper-triangle policy, while node4 must remain an isolated observation.
    row = [0, 0, 1, 1, 1, 2, 2, 3, 0, 3, 1]
    col = [1, 1, 0, 0, 2, 1, 3, 2, 3, 0, 1]
    data = np.array([.25, .5, .25, .5, 2.5, 2.5, .125, .125, 4.25, 4.25, 99.], dtype=dtype)
    graph = sparse.coo_matrix((data, (row, col)), shape=(5, 5))
    path = tmp_path / "nonuniform.npz"
    sparse.save_npz(path, graph)
    parallel_leiden._initialize(str(path), 1)
    native = parallel_leiden._GRAPH
    assert native.numberOfNodes() == 5
    assert native.numberOfEdges() == 4
    assert native.weight(0, 1) == .75
    assert native.weight(1, 2) == 2.5
    assert native.weight(2, 3) == .125
    assert native.weight(0, 3) == 4.25
    assert native.totalEdgeWeight() == 7.625
    np.testing.assert_array_equal([native.weightedDegree(i) for i in range(5)], [5., 3.25, 2.625, 4.375, 0.])
    assert not native.hasEdge(1, 1)


def test_legacy_builder_checkpoint_is_rejected_without_modification(tmp_path):
    import json

    pytest.importorskip("networkit")
    graph = sparse.csr_matrix([[0., .25, 0.], [.25, 0., 2.5], [0., 2.5, 0.]])
    config = {"checkpoint_dir": str(tmp_path), "workers": 1, "threads_per_worker": 1, "iterations": 2}
    list(run_jobs(graph, [(.5, 17)], config))
    metadata = tmp_path / "r0.5_s17.json"
    record = json.loads(metadata.read_text())
    assert record.pop("graph_build_policy") == parallel_leiden.GRAPH_BUILD_POLICY
    metadata.write_text(json.dumps(record))
    labels = tmp_path / "r0.5_s17.npz"
    metadata_before, labels_before = metadata.read_bytes(), labels.read_bytes()
    with pytest.raises(ValueError, match="Checkpoint settings changed"):
        list(run_jobs(graph, [(.5, 17)], config))
    assert metadata.read_bytes() == metadata_before
    assert labels.read_bytes() == labels_before


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
