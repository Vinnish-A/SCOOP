from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from .benchmark_manifest import FastCNMFBenchmarkManifest


@dataclass(frozen=True)
class ExecutionStage:
    """A concrete stage in an independent benchmark execution plan."""

    stage_id: str
    dataset_id: str
    lane_id: str
    kind: str
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]
    cold_start_required: bool
    may_reuse_cross_lane_artifacts: bool


@dataclass(frozen=True)
class ExecutionPlan:
    """Serializable execution plan derived from a benchmark manifest."""

    schema_version: str
    run_id: str
    stages: tuple[ExecutionStage, ...]

    def to_json(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2), encoding="utf-8")

    @classmethod
    def from_json(cls, path: Path) -> "ExecutionPlan":
        payload = json.loads(path.read_text(encoding="utf-8"))
        return cls(
            schema_version=payload["schema_version"],
            run_id=payload["run_id"],
            stages=tuple(ExecutionStage(**stage) for stage in payload["stages"]),
        )


def build_execution_plan(manifest: FastCNMFBenchmarkManifest) -> ExecutionPlan:
    """Expand a benchmark manifest into dataset/lane/stage tasks.

    This is the boundary between old wrapper scripts and the independent
    framework. The stages intentionally describe artifacts in FastCNMF-owned
    output roots so the candidate lane cannot accidentally consume cNMF
    reference intermediates.
    """

    stages: list[ExecutionStage] = []
    for dataset in manifest.datasets:
        for lane in manifest.lanes:
            base = Path(manifest.output_root) / dataset.dataset_id / lane.lane_id
            raw_input = dataset.path
            core = base / "preprocess/core_arrays.npz"
            preprocess_prefix = base / "preprocess/cnmf_input"
            corrected = Path(str(preprocess_prefix) + ".Corrected.HVG.Varnorm.h5ad")
            tp10k = Path(str(preprocess_prefix) + ".TP10K.h5ad")
            hvgs = Path(str(preprocess_prefix) + ".Corrected.HVGs.txt")
            task_manifest = base / "nmf/task_manifest.json"
            nmf_output = base / "nmf/replicates"
            consensus = base / "consensus"
            report = base / "reports/benchmark_summary.json"

            no_cross_lane_reuse = not lane.allow_existing_intermediates
            stages.extend(
                (
                    ExecutionStage(
                        stage_id=f"{dataset.dataset_id}:{lane.lane_id}:preprocess",
                        dataset_id=dataset.dataset_id,
                        lane_id=lane.lane_id,
                        kind="preprocess_cold_start",
                        inputs=(raw_input,),
                        outputs=(str(core), str(corrected), str(tp10k), str(hvgs)),
                        cold_start_required=lane.cold_start,
                        may_reuse_cross_lane_artifacts=not no_cross_lane_reuse,
                    ),
                    ExecutionStage(
                        stage_id=f"{dataset.dataset_id}:{lane.lane_id}:plan_nmf",
                        dataset_id=dataset.dataset_id,
                        lane_id=lane.lane_id,
                        kind="plan_nmf_tasks",
                        inputs=(str(corrected), str(tp10k), str(hvgs)),
                        outputs=(str(task_manifest),),
                        cold_start_required=lane.cold_start,
                        may_reuse_cross_lane_artifacts=not no_cross_lane_reuse,
                    ),
                    ExecutionStage(
                        stage_id=f"{dataset.dataset_id}:{lane.lane_id}:factorize",
                        dataset_id=dataset.dataset_id,
                        lane_id=lane.lane_id,
                        kind="factorize_replicates",
                        inputs=(str(corrected), str(task_manifest)),
                        outputs=(str(nmf_output),),
                        cold_start_required=lane.cold_start,
                        may_reuse_cross_lane_artifacts=not no_cross_lane_reuse,
                    ),
                    ExecutionStage(
                        stage_id=f"{dataset.dataset_id}:{lane.lane_id}:consensus",
                        dataset_id=dataset.dataset_id,
                        lane_id=lane.lane_id,
                        kind="consensus",
                        inputs=(str(nmf_output), str(tp10k)),
                        outputs=(str(consensus),),
                        cold_start_required=lane.cold_start,
                        may_reuse_cross_lane_artifacts=not no_cross_lane_reuse,
                    ),
                    ExecutionStage(
                        stage_id=f"{dataset.dataset_id}:{lane.lane_id}:report",
                        dataset_id=dataset.dataset_id,
                        lane_id=lane.lane_id,
                        kind="resource_and_quality_report",
                        inputs=(str(consensus),),
                        outputs=(str(report),),
                        cold_start_required=lane.cold_start,
                        may_reuse_cross_lane_artifacts=not no_cross_lane_reuse,
                    ),
                )
            )

    return ExecutionPlan(
        schema_version="fastcnmf.execution_plan.v1",
        run_id=manifest.run_id,
        stages=tuple(stages),
    )
