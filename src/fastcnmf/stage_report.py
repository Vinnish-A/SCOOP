from __future__ import annotations

import json
import re
from pathlib import Path


def parse_elapsed(value: str) -> float:
    parts = value.strip().split(":")
    if len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    return float(value)


def parse_time_log(path: Path) -> dict:
    result = {"path": str(path)}
    if not path.exists():
        result["missing"] = True
        return result
    text = path.read_text(errors="replace")
    signal_match = re.search(r"Command terminated by signal (\d+)", text)
    if signal_match:
        result["terminated_by_signal"] = int(signal_match.group(1))
    fields = {
        "Elapsed (wall clock) time (h:mm:ss or m:ss):": ("elapsed_seconds", parse_elapsed),
        "Maximum resident set size (kbytes):": ("max_rss_mb", lambda x: round(int(x) / 1024, 1)),
        "Percent of CPU this job got:": ("cpu_percent", lambda x: float(x.rstrip("%"))),
        "User time (seconds):": ("user_seconds", float),
        "System time (seconds):": ("system_seconds", float),
        "Exit status:": ("exit_status", int),
    }
    for line in text.splitlines():
        line = line.strip()
        for prefix, (key, fn) in fields.items():
            if line.startswith(prefix):
                result[key] = fn(line.removeprefix(prefix).strip())
                break
    return result


def directory_size_bytes(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file())


def write_stage_report(
    *,
    dataset_id: str,
    lane_id: str,
    lane_root: Path,
    logs_dir: Path,
    output_json: Path,
) -> dict:
    prefix = f"{dataset_id}__{lane_id}__"
    stages = {}
    for log in sorted(logs_dir.glob(f"{prefix}*.time.log")):
        stage_name = log.name.removeprefix(prefix).removesuffix(".time.log")
        stages[stage_name] = parse_time_log(log)

    report = {
        "dataset_id": dataset_id,
        "lane_id": lane_id,
        "lane_root": str(lane_root),
        "output_size_mb": round(directory_size_bytes(lane_root) / 1024**2, 2),
        "stages": stages,
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
