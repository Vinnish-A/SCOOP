from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping
import yaml


def read_yaml(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    with path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping: {path}")
    return data


def deep_get(config: Mapping[str, Any], dotted: str, default: Any = None) -> Any:
    cur: Any = config
    for key in dotted.split("."):
        if not isinstance(cur, Mapping) or key not in cur:
            return default
        cur = cur[key]
    return cur


def resolve_run_root(config_path: str | Path, config: Mapping[str, Any]) -> Path:
    """Return the run root from a config file path.

    The expected layout is runs/<run_id>/config/run.yaml. If a different
    layout is used, this function falls back to the parent directory of
    the config file.
    """
    config_path = Path(config_path).resolve()
    if config_path.parent.name == "config":
        return config_path.parent.parent
    run_id = deep_get(config, "run.run_id", "run")
    return Path("runs") / str(run_id)


def project_root_from_run_root(run_root: Path) -> Path:
    if run_root.parent.name == "runs":
        return run_root.parent.parent
    return Path.cwd()
