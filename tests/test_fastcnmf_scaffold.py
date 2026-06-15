from pathlib import Path

from fastcnmf.config import NMFConfig
from fastcnmf.fast_prepare import run_fast_prepare
from fastcnmf.resources import choose_replicate_batch_size, estimate_nmf_memory
from fastcnmf.tasks import TaskManifest, build_nmf_tasks
from scsp_agent_sop.programs import FASTCNMF_DEFAULT_MAX_ITER, FASTCNMF_DEFAULT_SEEDS


def test_build_nmf_tasks_roundtrip(tmp_path: Path) -> None:
    manifest = build_nmf_tasks("run", k_values=(6, 8), n_iter=3, seed=11)
    assert len(manifest.tasks) == 6
    assert manifest.tasks[0].task_id == "k6_iter0"

    path = tmp_path / "manifest.json"
    manifest.to_json(path)
    loaded = TaskManifest.from_json(path)
    assert loaded == manifest


def test_fastcnmf_defaults_are_scoop_preferred() -> None:
    assert FASTCNMF_DEFAULT_MAX_ITER == 50
    assert FASTCNMF_DEFAULT_SEEDS == tuple(range(20))
    assert NMFConfig().max_nmf_iter == 50
    assert run_fast_prepare.__kwdefaults__["max_nmf_iter"] == 50


def test_memory_estimate_and_batch_choice() -> None:
    estimate = estimate_nmf_memory(
        observations=10_000,
        genes=3_000,
        components=8,
        dtype="float32",
        replicate_batch_size=1,
    )
    assert estimate.estimated_bytes > 0
    assert estimate.estimated_gib > 0

    batch = choose_replicate_batch_size(
        observations=10_000,
        genes=3_000,
        components=8,
        available_bytes=8 * 1024**3,
        dtype="float32",
        max_batch_size=4,
    )
    assert 1 <= batch <= 4
