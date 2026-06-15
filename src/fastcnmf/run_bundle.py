from __future__ import annotations

import json
import shlex
from dataclasses import asdict, dataclass
from pathlib import Path

from .benchmark_manifest import FastCNMFBenchmarkManifest
from .runner import ExecutionPlan, ExecutionStage


@dataclass(frozen=True)
class StageScript:
    stage_id: str
    path: str
    implemented: bool
    cuda_visible_devices: str | None = None


@dataclass(frozen=True)
class RunBundle:
    """Collection of generated executable stage scripts."""

    schema_version: str
    run_id: str
    scripts: tuple[StageScript, ...]

    def to_json(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2), encoding="utf-8")


def _q(value: str | Path) -> str:
    return shlex.quote(str(value))


def _stage_script_name(stage: ExecutionStage) -> str:
    return stage.stage_id.replace(":", "__") + ".sh"


def _prefix_from_corrected_h5ad(path: str) -> Path:
    suffix = ".Corrected.HVG.Varnorm.h5ad"
    if not path.endswith(suffix):
        raise ValueError(f"corrected H5AD output must end with {suffix}: {path}")
    return Path(path[: -len(suffix)])


def _run_name(stage: ExecutionStage) -> str:
    return f"{stage.dataset_id}_{stage.lane_id}"


def _cnmf_binary(python: str) -> str:
    path = Path(python)
    if path.name.startswith("python"):
        return str(path.parent / "cnmf")
    return "cnmf"


def _nmf_root_from_task_manifest(path: str) -> Path:
    return Path(path).parent


def _cnmf_output_dir(stage: ExecutionStage) -> Path:
    if stage.kind == "plan_nmf_tasks":
        return _nmf_root_from_task_manifest(stage.outputs[0])
    if stage.kind == "factorize_replicates":
        return Path(stage.inputs[1]).parent
    if stage.kind == "consensus":
        return Path(stage.inputs[0]).parent
    if stage.kind == "resource_and_quality_report":
        return Path(stage.inputs[0]).parent.parent / "nmf"
    raise ValueError(f"cannot infer cNMF output dir for stage {stage.stage_id}")


def _profile_value(params: dict, key: str, profile: str):
    defaults = {
        "n_iter": 50,
        "max_nmf_iter": 50,
    }
    profiled = f"{key}_{profile}"
    if profiled in params:
        return params[profiled]
    return params.get(key, defaults[key])


def _preprocess_command(stage: ExecutionStage, manifest: FastCNMFBenchmarkManifest, python: str) -> tuple[str, bool]:
    params = manifest.parameters
    output_prefix = _prefix_from_corrected_h5ad(stage.outputs[1])
    if stage.lane_id == "cnmf_optimized":
        return (
            " ".join(
                [
                    _q(python),
                    "-m fastcnmf cnmf-preprocess",
                    "--input-h5ad",
                    _q(stage.inputs[0]),
                    "--output-prefix",
                    _q(output_prefix),
                    "--sample-key",
                    _q(params.get("harmony_batch_key", "sample_id")),
                    "--n-top-genes",
                    str(params.get("n_top_genes", 3000)),
                    "--theta",
                    "1",
                    "--max-iter-harmony",
                    str(params.get("max_iter_harmony", 20)),
                    "--seed",
                    str(params.get("seed", 20260614)),
                ]
            ),
            True,
        )
    return (
        " ".join(
            [
                _q(python),
                "-m fastcnmf fast-preprocess",
                "--input-h5ad",
                _q(stage.inputs[0]),
                "--output-prefix",
                _q(output_prefix),
                "--sample-key",
                _q(params.get("harmony_batch_key", "sample_id")),
                "--n-top-genes",
                str(params.get("n_top_genes", 3000)),
                "--theta",
                "1",
                "--max-iter-harmony",
                str(params.get("max_iter_harmony", 20)),
                "--seed",
                str(params.get("seed", 20260614)),
            ]
        ),
        True,
    )


def _plan_nmf_command(
    stage: ExecutionStage,
    manifest: FastCNMFBenchmarkManifest,
    python: str,
    *,
    profile: str,
) -> tuple[str, bool]:
    params = manifest.parameters
    cnmf = _cnmf_binary(python)
    corrected, tp10k, hvgs = stage.inputs
    output_dir = _cnmf_output_dir(stage)
    n_iter = int(_profile_value(params, "n_iter", profile))
    max_nmf_iter = int(_profile_value(params, "max_nmf_iter", profile))
    k_values = [str(k) for k in params.get("k_values", [6, 8, 10, 12])]
    run_name = _run_name(stage)
    if stage.lane_id == manifest.fairness.candidate_lane:
        return (
            " && ".join(
                [
                    " ".join(
                        [
                            _q(python),
                            "-m fastcnmf plan-tasks",
                            "--run-name",
                            _q(run_name),
                            "-k",
                            *k_values,
                            "--n-iter",
                            str(n_iter),
                            "--seed",
                            str(params.get("seed", 20260614)),
                            "--output",
                            _q(stage.outputs[0]),
                        ]
                    ),
                    " ".join(
                        [
                            _q(python),
                            "-m fastcnmf fast-prepare",
                            "--corrected-h5ad",
                            _q(corrected),
                            "--tpm-h5ad",
                            _q(tp10k),
                            "--hvgs-txt",
                            _q(hvgs),
                            "--output-dir",
                            _q(output_dir),
                            "--run-name",
                            _q(run_name),
                            "-k",
                            *k_values,
                            "--n-iter",
                            str(n_iter),
                            "--seed",
                            str(params.get("seed", 20260614)),
                            "--max-nmf-iter",
                            str(max_nmf_iter),
                        ]
                    ),
                ]
            ),
            True,
        )

    return (
        " && ".join(
            [
                " ".join(
                    [
                        _q(python),
                        "-m fastcnmf plan-tasks",
                        "--run-name",
                        _q(run_name),
                        "-k",
                        *k_values,
                        "--n-iter",
                        str(n_iter),
                        "--seed",
                        str(params.get("seed", 20260614)),
                        "--output",
                        _q(stage.outputs[0]),
                    ]
                ),
                " ".join(
                    [
                        _q(cnmf),
                        "prepare",
                        "--output-dir",
                        _q(output_dir),
                        "--name",
                        _q(run_name),
                        "--counts",
                        _q(corrected),
                        "--tpm",
                        _q(tp10k),
                        "--genes-file",
                        _q(hvgs),
                        "-k",
                        *k_values,
                        "--n-iter",
                        str(n_iter),
                        "--seed",
                        str(params.get("seed", 20260614)),
                        "--numgenes",
                        str(params.get("n_top_genes", 3000)),
                        "--max-nmf-iter",
                        str(max_nmf_iter),
                    ]
                ),
            ]
        ),
        True,
    )


def _factorize_command(stage: ExecutionStage, manifest: FastCNMFBenchmarkManifest, python: str) -> tuple[str, bool]:
    params = manifest.parameters
    output_dir = _cnmf_output_dir(stage)
    run_name = _run_name(stage)
    workers = int(params.get("workers", 4))
    if stage.lane_id == manifest.fairness.candidate_lane:
        return (
            " ".join(
                [
                    _q(python),
                    "-m fastcnmf fast-factorize",
                    "--output-dir",
                    _q(output_dir),
                    "--run-name",
                    _q(run_name),
                    "--workers",
                    str(workers),
                ]
            ),
            True,
        )

    cnmf = _cnmf_binary(python)
    worker_log_dir = output_dir / "worker_logs"
    snippet = f"""
mkdir -p {_q(worker_log_dir)}
status=0
for worker_i in $(seq 0 {workers - 1}); do
  (
    {_q(cnmf)} factorize --output-dir {_q(output_dir)} --name {_q(run_name)} --total-workers {workers} --worker-index "$worker_i" \
      > {_q(worker_log_dir)}/factorize_worker_${{worker_i}}.stdout.log \
      2> {_q(worker_log_dir)}/factorize_worker_${{worker_i}}.stderr.log
  ) &
done
for job in $(jobs -p); do
  if ! wait "$job"; then
    status=1
  fi
done
exit "$status"
""".strip()
    return snippet, True


def _consensus_command(stage: ExecutionStage, manifest: FastCNMFBenchmarkManifest, python: str) -> tuple[str, bool]:
    params = manifest.parameters
    output_dir = _cnmf_output_dir(stage)
    run_name = _run_name(stage)
    if stage.lane_id == manifest.fairness.candidate_lane:
        workers = int(params.get("workers", 4))
        k_values = [str(k) for k in params.get("k_values", [6, 8, 10, 12])]
        return (
            " ".join(
                [
                    _q(python),
                    "-m fastcnmf fast-consensus",
                    "--output-dir",
                    _q(output_dir),
                    "--run-name",
                    _q(run_name),
                    "--workers",
                    str(workers),
                    "-k",
                    *k_values,
                ]
            ),
            True,
        )

    cnmf = _cnmf_binary(python)
    return (
        " && ".join(
            [
                f"{_q(cnmf)} combine --output-dir {_q(output_dir)} --name {_q(run_name)}",
                f"{_q(cnmf)} consensus --output-dir {_q(output_dir)} --name {_q(run_name)}",
                f"{_q(cnmf)} k_selection_plot --output-dir {_q(output_dir)} --name {_q(run_name)}",
            ]
        ),
        True,
    )


def _report_command(stage: ExecutionStage, manifest: FastCNMFBenchmarkManifest, python: str, logs_dir: Path) -> tuple[str, bool]:
    lane_root = Path(manifest.output_root) / stage.dataset_id / stage.lane_id
    return (
        " ".join(
            [
                _q(python),
                "-m fastcnmf stage-report",
                "--dataset-id",
                _q(stage.dataset_id),
                "--lane-id",
                _q(stage.lane_id),
                "--lane-root",
                _q(lane_root),
                "--logs-dir",
                _q(logs_dir),
                "--output-json",
                _q(stage.outputs[0]),
            ]
        ),
        True,
    )


def _placeholder_command(stage: ExecutionStage) -> tuple[str, bool]:
    msg = (
        f"Stage executor for {stage.kind!r} / lane {stage.lane_id!r} is not implemented yet. "
        "The script is generated to preserve the benchmark artifact boundary."
    )
    return f"printf '%s\\n' {_q(msg)}; exit 78", False


def stage_command(
    stage: ExecutionStage,
    manifest: FastCNMFBenchmarkManifest,
    *,
    reference_python: str,
    candidate_python: str,
    profile: str,
    logs_dir: Path,
) -> tuple[str, bool]:
    python = candidate_python if stage.lane_id == manifest.fairness.candidate_lane else reference_python
    if stage.kind == "preprocess_cold_start":
        return _preprocess_command(stage, manifest, python)
    if stage.kind == "plan_nmf_tasks":
        return _plan_nmf_command(stage, manifest, python, profile=profile)
    if stage.kind == "factorize_replicates":
        return _factorize_command(stage, manifest, python)
    if stage.kind == "consensus":
        return _consensus_command(stage, manifest, python)
    if stage.kind == "resource_and_quality_report":
        return _report_command(stage, manifest, python, logs_dir)
    return _placeholder_command(stage)


def write_run_bundle(
    *,
    manifest: FastCNMFBenchmarkManifest,
    plan: ExecutionPlan,
    output_dir: Path,
    reference_python: str = "python",
    candidate_python: str = "python",
    reference_cuda_visible_devices: str | None = None,
    candidate_cuda_visible_devices: str | None = None,
    profile: str = "production",
    blas_threads: int = 1,
) -> RunBundle:
    scripts_dir = output_dir / "scripts"
    logs_dir = output_dir / "logs"
    scripts_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    scripts: list[StageScript] = []
    for stage in plan.stages:
        cuda_visible_devices = (
            candidate_cuda_visible_devices
            if stage.lane_id == manifest.fairness.candidate_lane
            else reference_cuda_visible_devices
        )
        command, implemented = stage_command(
            stage,
            manifest,
            reference_python=reference_python,
            candidate_python=candidate_python,
            profile=profile,
            logs_dir=logs_dir,
        )
        script_path = scripts_dir / _stage_script_name(stage)
        log_base = logs_dir / stage.stage_id.replace(":", "__")
        for output in stage.outputs:
            Path(output).parent.mkdir(parents=True, exist_ok=True)
        script_path.write_text(
            "\n".join(
                [
                    "#!/usr/bin/env bash",
                    "set -euo pipefail",
                    f"export OMP_NUM_THREADS={blas_threads}",
                    f"export MKL_NUM_THREADS={blas_threads}",
                    f"export OPENBLAS_NUM_THREADS={blas_threads}",
                    f"export NUMEXPR_NUM_THREADS={blas_threads}",
                    'export PYTHONPATH="${PYTHONPATH:-}:src"',
                    (
                        f"export CUDA_VISIBLE_DEVICES={shlex.quote(cuda_visible_devices)}"
                        if cuda_visible_devices is not None
                        else "# CUDA_VISIBLE_DEVICES inherited from parent environment"
                    ),
                    f"mkdir -p {_q(log_base.parent)}",
                    f"echo {shlex.quote(stage.stage_id)} > {_q(log_base)}.stage_id",
                    f"/usr/bin/time -v -o {_q(log_base)}.time.log bash -lc {_q(command)} "
                    f"> {_q(log_base)}.stdout.log 2> {_q(log_base)}.stderr.log",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        script_path.chmod(0o755)
        scripts.append(
            StageScript(
                stage_id=stage.stage_id,
                path=str(script_path),
                implemented=implemented,
                cuda_visible_devices=cuda_visible_devices,
            )
        )

    return RunBundle(
        schema_version="fastcnmf.run_bundle.v1",
        run_id=plan.run_id,
        scripts=tuple(scripts),
    )
