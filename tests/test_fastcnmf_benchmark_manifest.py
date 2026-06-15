from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from fastcnmf.benchmark_manifest import (
    BenchmarkLane,
    FastCNMFBenchmarkManifest,
    FairnessPolicy,
    inspect_h5ad_dataset,
)
from fastcnmf.runner import build_execution_plan
from fastcnmf.run_bundle import write_run_bundle


def _write_h5ad(path: Path) -> None:
    obs = pd.DataFrame(
        {"sample_id": ["s1", "s1", "s2", "s2", "s2", "s3"]},
        index=[f"cell{i}" for i in range(6)],
    )
    var = pd.DataFrame(index=[f"gene{i}" for i in range(4)])
    ad.AnnData(np.ones((6, 4)), obs=obs, var=var).write_h5ad(path)


def test_inspect_h5ad_dataset_records_sample_scale(tmp_path: Path) -> None:
    h5ad = tmp_path / "dataset.h5ad"
    _write_h5ad(h5ad)

    spec = inspect_h5ad_dataset(h5ad, dataset_id="s2_test", tier="S2", modality="single_cell")

    assert spec.n_obs == 6
    assert spec.n_vars == 4
    assert spec.sample_key == "sample_id"
    assert spec.n_samples == 3
    assert spec.min_cells_per_sample == 1
    assert spec.max_cells_per_sample == 3


def test_execution_plan_keeps_candidate_lane_cold_start_isolated(tmp_path: Path) -> None:
    # Construct directly here so this test focuses on plan semantics.
    from fastcnmf.benchmark_manifest import DatasetSpec

    manifest = FastCNMFBenchmarkManifest(
        schema_version="fastcnmf.benchmark_manifest.v1",
        run_id="test",
        generated_from=str(tmp_path),
        datasets=(
            DatasetSpec(
                dataset_id="s2_test",
                tier="S2",
                modality="single_cell",
                path=str(tmp_path / "dataset.h5ad"),
                source_kind="h5ad",
                n_obs=72_000,
                n_vars=30_000,
            ),
        ),
        lanes=(
            BenchmarkLane(
                lane_id="cnmf_optimized",
                engine="cnmf",
                cold_start=True,
                allow_existing_intermediates=False,
                allowed_hardware_acceleration=("multi_process_cpu",),
                description="baseline",
            ),
            BenchmarkLane(
                lane_id="fastcnmf_independent",
                engine="fastcnmf",
                cold_start=True,
                allow_existing_intermediates=False,
                allowed_hardware_acceleration=("multi_process_cpu", "gpu_if_enabled"),
                description="candidate",
            ),
        ),
        fairness=FairnessPolicy(target_speedup=3.0),
    )

    plan = build_execution_plan(manifest)

    assert len(plan.stages) == 10
    candidate = [stage for stage in plan.stages if stage.lane_id == "fastcnmf_independent"]
    assert [stage.kind for stage in candidate] == [
        "preprocess_cold_start",
        "plan_nmf_tasks",
        "factorize_replicates",
        "consensus",
        "resource_and_quality_report",
    ]
    assert all(stage.cold_start_required for stage in candidate)
    assert not any(stage.may_reuse_cross_lane_artifacts for stage in candidate)


def test_run_bundle_writes_reference_preprocess_and_candidate_placeholder(tmp_path: Path) -> None:
    from fastcnmf.benchmark_manifest import DatasetSpec

    manifest = FastCNMFBenchmarkManifest(
        schema_version="fastcnmf.benchmark_manifest.v1",
        run_id="test",
        generated_from=str(tmp_path),
        datasets=(
            DatasetSpec(
                dataset_id="s2_test",
                tier="S2",
                modality="single_cell",
                path=str(tmp_path / "dataset.h5ad"),
                source_kind="h5ad",
                n_obs=72_000,
                n_vars=30_000,
            ),
        ),
        lanes=(
            BenchmarkLane(
                lane_id="cnmf_optimized",
                engine="cnmf",
                cold_start=True,
                allow_existing_intermediates=False,
                allowed_hardware_acceleration=("multi_process_cpu",),
                description="baseline",
            ),
            BenchmarkLane(
                lane_id="fastcnmf_independent",
                engine="fastcnmf",
                cold_start=True,
                allow_existing_intermediates=False,
                allowed_hardware_acceleration=("multi_process_cpu",),
                description="candidate",
            ),
        ),
        fairness=FairnessPolicy(target_speedup=3.0),
        output_root=str(tmp_path / "benchmark"),
    )
    plan = build_execution_plan(manifest)

    bundle = write_run_bundle(
        manifest=manifest,
        plan=plan,
        output_dir=tmp_path / "bundle",
        reference_python="python-ref",
        candidate_python="python-fast",
        reference_cuda_visible_devices="",
        candidate_cuda_visible_devices="0",
    )

    ref_preprocess = next(script for script in bundle.scripts if script.stage_id == "s2_test:cnmf_optimized:preprocess")
    cand_preprocess = next(script for script in bundle.scripts if script.stage_id == "s2_test:fastcnmf_independent:preprocess")
    assert ref_preprocess.implemented
    assert cand_preprocess.implemented
    ref_text = Path(ref_preprocess.path).read_text(encoding="utf-8")
    cand_text = Path(cand_preprocess.path).read_text(encoding="utf-8")
    assert "cnmf-preprocess" in ref_text
    assert "python-ref" in ref_text
    assert "CUDA_VISIBLE_DEVICES=" in ref_text
    assert "PYTHONPATH" in ref_text
    assert "fast-preprocess" in cand_text
    assert "python-fast" in cand_text
    assert "CUDA_VISIBLE_DEVICES=0" in cand_text

    cand_plan = next(script for script in bundle.scripts if script.stage_id == "s2_test:fastcnmf_independent:plan_nmf")
    cand_factorize = next(script for script in bundle.scripts if script.stage_id == "s2_test:fastcnmf_independent:factorize")
    cand_consensus = next(script for script in bundle.scripts if script.stage_id == "s2_test:fastcnmf_independent:consensus")
    cand_plan_text = Path(cand_plan.path).read_text(encoding="utf-8")
    assert "fast-prepare" in cand_plan_text
    assert "cnmf prepare" not in cand_plan_text
    assert "fast-factorize" in Path(cand_factorize.path).read_text(encoding="utf-8")
    assert "fast-consensus" in Path(cand_consensus.path).read_text(encoding="utf-8")
