#!/usr/bin/env python
"""Summarize the cNMF Harmony spatial benchmark."""
from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TMP = ROOT / "tmp/cnmf_spatial_harmony_benchmark"
LOGS = TMP / "logs"
OUT_JSON = TMP / "benchmark_summary.json"
OUT_MD = TMP / "benchmark_report.md"


def parse_elapsed(value: str) -> float:
    parts = value.strip().split(":")
    if len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    return float(value)


def parse_time_log(path: Path) -> dict:
    data: dict = {"path": str(path)}
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
    for line in path.read_text(errors="replace").splitlines():
        line = line.strip()
        for prefix, (key, converter) in fields.items():
            if line.startswith(prefix):
                data[key] = converter(line.removeprefix(prefix).strip())
                break
    return data


def du_mb(path: Path) -> float:
    return round(sum(p.stat().st_size for p in path.rglob("*") if p.is_file()) / 1024 / 1024, 1)


def collect() -> dict:
    manifest = json.loads((TMP / "input/gbm_lowres_visium_3samples_harmony_manifest.json").read_text())
    summary = {
        "input": manifest,
        "prepare_input": parse_time_log(TMP / "prepare_input.time.log"),
        "environment": {
            "venv": str(ROOT / ".venv-cnmf"),
            "cnmf_version": "1.7.1",
            "harmonypy_version": "0.2.0",
            "python": str(ROOT / ".venv-cnmf/bin/python"),
            "thread_limits": {
                "OMP_NUM_THREADS": 1,
                "MKL_NUM_THREADS": 1,
                "OPENBLAS_NUM_THREADS": 1,
                "NUMEXPR_NUM_THREADS": 1,
            },
        },
        "parameters": {
            "harmony_vars": ["sample_id"],
            "n_top_rna_genes": 3000,
            "k_values": [6, 8],
            "n_iter": 8,
            "total_nmf_runs": 16,
            "max_nmf_iter": 200,
            "parallel_workers": 4,
        },
        "runs": {},
    }

    for mode in ["serial", "parallel"]:
        run = {"output_dir": str(TMP / mode), "output_size_mb": du_mb(TMP / mode), "steps": {}}
        for step in ["prepare", "factorize", "combine", "consensus", "k_selection_plot"]:
            time_log = LOGS / mode / f"{step}.time.log"
            if time_log.exists():
                run["steps"][step] = parse_time_log(time_log)
        if mode == "parallel":
            text = (LOGS / mode / "factorize_parallel.cmd.log").read_text()
            elapsed = int(re.search(r"elapsed_seconds=(\d+)", text).group(1))
            worker_logs = {}
            for worker_i in range(4):
                worker_logs[str(worker_i)] = parse_time_log(LOGS / mode / f"factorize_worker_{worker_i}.time.log")
            run["steps"]["factorize_parallel_group"] = {
                "workers": 4,
                "elapsed_seconds": elapsed,
                "worker_logs": worker_logs,
            }
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


def write_md(summary: dict) -> None:
    def row(mode: str, step: str, label: str | None = None) -> str:
        data = summary["runs"][mode]["steps"][step]
        return (
            f"| {mode} | {label or step} | {data.get('elapsed_seconds', '')} | "
            f"{data.get('cpu_percent', '')} | {data.get('max_rss_mb', '')} | "
            f"{data.get('fs_inputs', '')} | {data.get('fs_outputs', '')} |"
        )

    worker_rows = "\n".join(
        f"| {idx} | {d['elapsed_seconds']} | {d['cpu_percent']} | {d['max_rss_mb']} | {d['exit_status']} |"
        for idx, d in summary["runs"]["parallel"]["steps"]["factorize_parallel_group"]["worker_logs"].items()
    )
    inp = summary["input"]
    cmp = summary["comparison"]
    md = f"""# Spatial cNMF Harmony Benchmark

## Input and Harmony Preprocessing

- Merged raw H5AD: `{inp['merged_raw_h5ad']}`
- Harmony-corrected cNMF input: `{inp['corrected_h5ad']}`
- TP10K input passed to cNMF: `{inp['tp10k_h5ad']}`
- HVG list passed to cNMF: `{inp['hvg_txt']}`
- Batch variable removed by Harmony: `{inp['batch_key_removed']}`
- Method: `{inp['method']}`
- Spots: `{inp['n_obs']}`
- Genes: `{inp['n_vars']}`
- Samples: `{', '.join(inp['samples'])}`
- Harmony preprocessing time: `{summary['prepare_input']['elapsed_seconds']} s`
- Harmony preprocessing max RSS: `{summary['prepare_input']['max_rss_mb']} MB`

## Environment

- UV venv: `{summary['environment']['venv']}`
- cNMF: `{summary['environment']['cnmf_version']}`
- harmonypy: `{summary['environment']['harmonypy_version']}`
- Thread limits: `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`,
  `OPENBLAS_NUM_THREADS=1`, `NUMEXPR_NUM_THREADS=1`

## cNMF Parameters

- k values: `{summary['parameters']['k_values']}`
- NMF replicates per k: `{summary['parameters']['n_iter']}`
- Total NMF runs: `{summary['parameters']['total_nmf_runs']}`
- max NMF iterations: `{summary['parameters']['max_nmf_iter']}`
- parallel workers: `{summary['parameters']['parallel_workers']}`

## Runtime and Resource Usage

| mode | step | elapsed_seconds | cpu_percent | max_rss_mb | fs_inputs | fs_outputs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
{row('serial', 'prepare')}
{row('serial', 'factorize')}
{row('serial', 'combine')}
{row('serial', 'consensus')}
{row('serial', 'k_selection_plot')}
{row('parallel', 'prepare')}
| parallel | factorize_group_4_workers | {summary['runs']['parallel']['steps']['factorize_parallel_group']['elapsed_seconds']} |  |  |  |  |
{row('parallel', 'combine')}
{row('parallel', 'consensus')}
{row('parallel', 'k_selection_plot')}

Output sizes:

- serial: `{summary['runs']['serial']['output_size_mb']} MB`
- parallel: `{summary['runs']['parallel']['output_size_mb']} MB`

## Parallel Worker Logs

| worker_index | elapsed_seconds | cpu_percent | max_rss_mb | exit_status |
| ---: | ---: | ---: | ---: | ---: |
{worker_rows}

## Serial vs Parallel

- Serial factorize: `{cmp['serial_factorize_seconds']} s`
- Parallel factorize, 4 workers: `{cmp['parallel_factorize_seconds']} s`
- Speedup: `{cmp['factorize_speedup']}x`
- Time saved: `{cmp['factorize_time_saved_seconds']} s`

## Parallelization Method

The cNMF factorization step was parallelized using four independent
`cnmf factorize` processes with `--total-workers 4` and worker indices `0..3`.
The Harmony preprocessing, cNMF prepare, combine, consensus, and k-selection
steps were run once per output directory.

## Issues Encountered

- `harmonypy==2.0.0` is not compatible with `cnmf==1.7.1` preprocessing because
  cNMF expects `Harmony.Phi_moe`. The environment was switched to
  `harmonypy==0.2.0`, which provides `Phi_moe`.
- Installing `harmonypy<2` pulled large Torch/CUDA dependencies into the UV
  environment. This affects environment size, not the recorded cNMF runtime.
- cNMF CLI `prepare` does not expose `harmony_vars`; Harmony preprocessing must
  be run via `cnmf.preprocess.Preprocess` before the standard CLI workflow.
- The corrected input and TP10K/HVG files must be passed together:
  `--counts Corrected.HVG.Varnorm.h5ad --tpm TP10K.h5ad --genes-file Corrected.HVGs.txt`.
"""
    OUT_MD.write_text(md, encoding="utf-8")


def main() -> None:
    summary = collect()
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_md(summary)
    print(json.dumps(summary["comparison"], indent=2))
    print(f"Wrote {OUT_JSON}")
    print(f"Wrote {OUT_MD}")


if __name__ == "__main__":
    main()
