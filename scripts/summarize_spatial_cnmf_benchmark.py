#!/usr/bin/env python
"""Summarize spatial cNMF benchmark timing/resource logs."""
from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TMP = ROOT / "tmp/cnmf_spatial_benchmark"
LOGS = TMP / "logs"
OUT_JSON = TMP / "benchmark_summary.json"
OUT_MD = TMP / "benchmark_report.md"


def parse_elapsed(value: str) -> float:
    value = value.strip()
    parts = value.split(":")
    if len(parts) == 2:
        minutes, seconds = parts
        return int(minutes) * 60 + float(seconds)
    if len(parts) == 3:
        hours, minutes, seconds = parts
        return int(hours) * 3600 + int(minutes) * 60 + float(seconds)
    return float(value)


def parse_time_log(path: Path) -> dict[str, float | int | str]:
    data: dict[str, float | int | str] = {"path": str(path)}
    for line in path.read_text(errors="replace").splitlines():
        line = line.strip()
        fields = {
            "User time (seconds):": ("user_seconds", float),
            "System time (seconds):": ("system_seconds", float),
            "Percent of CPU this job got:": ("cpu_percent", lambda x: float(x.rstrip("%"))),
            "Elapsed (wall clock) time (h:mm:ss or m:ss):": ("elapsed_seconds", parse_elapsed),
            "Maximum resident set size (kbytes):": ("max_rss_mb", lambda x: round(int(x) / 1024, 1)),
            "File system inputs:": ("fs_inputs", int),
            "File system outputs:": ("fs_outputs", int),
            "Exit status:": ("exit_status", int),
        }
        for prefix, (key, converter) in fields.items():
            if line.startswith(prefix):
                data[key] = converter(line.removeprefix(prefix).strip())
                break
    return data


def du_mb(path: Path) -> float:
    total = 0
    for item in path.rglob("*"):
        if item.is_file():
            total += item.stat().st_size
    return round(total / 1024 / 1024, 1)


def collect() -> dict:
    summary: dict = {
        "environment": {
            "venv": str(ROOT / ".venv-cnmf"),
            "cnmf_version": "1.7.1",
            "python": str(ROOT / ".venv-cnmf/bin/python"),
            "thread_limits": {
                "OMP_NUM_THREADS": 1,
                "MKL_NUM_THREADS": 1,
                "OPENBLAS_NUM_THREADS": 1,
                "NUMEXPR_NUM_THREADS": 1,
            },
        },
        "input": json.loads((TMP / "input/gbm_lowres_visium_3samples_cnmf_input_manifest.json").read_text()),
        "parameters": {
            "k_values": [6, 8],
            "n_iter": 8,
            "total_nmf_runs": 16,
            "max_nmf_iter": 200,
            "seed": 20260614,
            "parallel_workers": 4,
        },
        "runs": {},
    }

    for mode in ["serial", "parallel"]:
        mode_logs = LOGS / mode
        run: dict = {
            "output_dir": str(TMP / mode),
            "output_size_mb": du_mb(TMP / mode),
            "steps": {},
        }
        for step in ["prepare", "factorize", "combine", "consensus", "k_selection_plot"]:
            time_log = mode_logs / f"{step}.time.log"
            if time_log.exists():
                run["steps"][step] = parse_time_log(time_log)
        if mode == "parallel":
            cmd_log = mode_logs / "factorize_parallel.cmd.log"
            text = cmd_log.read_text() if cmd_log.exists() else ""
            elapsed_match = re.search(r"elapsed_seconds=(\d+)", text)
            run["steps"]["factorize_parallel_group"] = {
                "workers": 4,
                "elapsed_seconds": int(elapsed_match.group(1)) if elapsed_match else None,
                "worker_logs": {},
            }
            for worker_i in range(4):
                worker_log = mode_logs / f"factorize_worker_{worker_i}.time.log"
                if worker_log.exists():
                    run["steps"]["factorize_parallel_group"]["worker_logs"][str(worker_i)] = parse_time_log(worker_log)
        summary["runs"][mode] = run

    serial_factorize = summary["runs"]["serial"]["steps"]["factorize"]["elapsed_seconds"]
    parallel_factorize = summary["runs"]["parallel"]["steps"]["factorize_parallel_group"]["elapsed_seconds"]
    summary["comparison"] = {
        "serial_factorize_seconds": serial_factorize,
        "parallel_factorize_seconds": parallel_factorize,
        "factorize_speedup": round(serial_factorize / parallel_factorize, 2),
        "factorize_time_saved_seconds": round(serial_factorize - parallel_factorize, 2),
    }
    return summary


def write_markdown(summary: dict) -> None:
    serial = summary["runs"]["serial"]
    parallel = summary["runs"]["parallel"]
    cmp = summary["comparison"]
    inp = summary["input"]

    def step_row(mode: str, step: str, label: str | None = None) -> str:
        data = summary["runs"][mode]["steps"][step]
        return (
            f"| {mode} | {label or step} | {data.get('elapsed_seconds', '')} | "
            f"{data.get('cpu_percent', '')} | {data.get('max_rss_mb', '')} | "
            f"{data.get('fs_inputs', '')} | {data.get('fs_outputs', '')} |"
        )

    worker_rows = []
    for worker_i, data in parallel["steps"]["factorize_parallel_group"]["worker_logs"].items():
        worker_rows.append(
            f"| {worker_i} | {data.get('elapsed_seconds')} | {data.get('cpu_percent')} | "
            f"{data.get('max_rss_mb')} | {data.get('exit_status')} |"
        )

    md = f"""# Spatial cNMF benchmark

## Environment

- UV venv: `{summary['environment']['venv']}`
- cNMF: `{summary['environment']['cnmf_version']}`
- Python: `{summary['environment']['python']}`
- Thread limits for both serial and parallel runs: `OMP_NUM_THREADS=1`,
  `MKL_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `NUMEXPR_NUM_THREADS=1`

## Input

- H5AD: `{inp['output_h5ad']}`
- Spots: `{inp['n_obs']}`
- Genes: `{inp['n_vars']}`
- Samples: `{', '.join(inp['samples'])}`
- Spots by sample: `{inp['spots_by_sample']}`

## cNMF Parameters

- k values: `{summary['parameters']['k_values']}`
- NMF replicates per k: `{summary['parameters']['n_iter']}`
- Total NMF runs: `{summary['parameters']['total_nmf_runs']}`
- max NMF iterations: `{summary['parameters']['max_nmf_iter']}`
- seed: `{summary['parameters']['seed']}`
- parallel workers: `{summary['parameters']['parallel_workers']}`

## Runtime and Resource Usage

| mode | step | elapsed_seconds | cpu_percent | max_rss_mb | fs_inputs | fs_outputs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
{step_row('serial', 'prepare')}
{step_row('serial', 'factorize')}
{step_row('serial', 'combine')}
{step_row('serial', 'consensus')}
{step_row('serial', 'k_selection_plot')}
{step_row('parallel', 'prepare')}
| parallel | factorize_group_4_workers | {parallel['steps']['factorize_parallel_group']['elapsed_seconds']} |  |  |  |  |
{step_row('parallel', 'combine')}
{step_row('parallel', 'consensus')}
{step_row('parallel', 'k_selection_plot')}

Output sizes:

- serial: `{serial['output_size_mb']} MB`
- parallel: `{parallel['output_size_mb']} MB`

## Parallel Factorize Worker Logs

| worker_index | elapsed_seconds | cpu_percent | max_rss_mb | exit_status |
| ---: | ---: | ---: | ---: | ---: |
{chr(10).join(worker_rows)}

## Serial vs Parallel

- Serial factorize: `{cmp['serial_factorize_seconds']} s`
- Parallel factorize, 4 workers: `{cmp['parallel_factorize_seconds']} s`
- Speedup: `{cmp['factorize_speedup']}x`
- Time saved: `{cmp['factorize_time_saved_seconds']} s`

## Parallelization Method

cNMF was parallelized by running four independent `cnmf factorize` processes
with `--total-workers 4` and `--worker-index 0..3`. `prepare`, `combine`,
`consensus`, and `k_selection_plot` were run once per output directory and were
not parallelized.

## Issues Encountered

- `uv venv` created a venv without an importable `pip` module. Installing with
  `uv pip install --python .venv-cnmf/bin/python ...` worked.
- Scanpy's `seurat_v3` HVG selection required `scikit-misc`; it was installed
  into the UV environment before preparing the H5AD input.
- To avoid nested BLAS oversubscription during the cNMF worker comparison,
  BLAS/OpenMP thread counts were pinned to one thread per worker.
- Parallel worker stdout/stderr had to be separated per worker to avoid
  interleaved logs.
"""
    OUT_MD.write_text(md, encoding="utf-8")


def main() -> None:
    summary = collect()
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_markdown(summary)
    print(json.dumps(summary["comparison"], indent=2))
    print(f"Wrote {OUT_JSON}")
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
