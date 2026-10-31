from __future__ import annotations

from dataclasses import dataclass
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
    return config_path.parent


def project_root_from_run_root(run_root: Path) -> Path:
    if run_root.parent.name == "runs":
        return run_root.parent.parent
    return Path.cwd()


OBS_KEY_ROLES = ("sample", "donor", "batch", "study")
# Pre-R10 config paths and built-in names, used when run.obs_keys leaves a role null.
_LEGACY_OBS_KEYS = {"sample": ("keys.sample", "sample_id"), "donor": ("keys.donor", "donor_id"),
                    "batch": ("keys.batch_candidates", None), "study": ("keys.study", "study_id")}


@dataclass(frozen=True)
class ObsKeys:
    """Typed obs columns (R10) and the config path each name came from."""

    sample: str | None
    donor: str | None
    batch: str | None
    study: str | None
    sources: Mapping[str, str]

    def record(self) -> dict[str, dict[str, str | None]]:
        return {role: {"key": getattr(self, role), "source": self.sources[role]} for role in OBS_KEY_ROLES}


def resolve_obs_keys(
    config: Mapping[str, Any],
    obs=None,
    *,
    overrides: Mapping[str, str] | None = None,
    required: tuple[str, ...] = (),
) -> ObsKeys:
    """Resolve ``run.obs_keys`` once, falling back to the legacy ``keys.*`` names.

    ``overrides`` maps a role to a stage-specific dotted config path that wins
    when set. With ``obs``, a declared key (override or ``run.obs_keys``) and
    every ``required`` role must be an obs column; a legacy or default name
    may be absent, as before. The batch fallback is the first
    ``keys.batch_candidates`` column present in ``obs``.
    """
    values: dict[str, str | None] = {}
    sources: dict[str, str] = {}
    declared = set()
    for role in OBS_KEY_ROLES:
        paths = [path for path in ((overrides or {}).get(role), f"run.obs_keys.{role}") if path]
        path = next((path for path in paths if deep_get(config, path) is not None), None)
        if path is not None:
            values[role], sources[role] = str(deep_get(config, path)), path
            declared.add(role)
            continue
        legacy, default = _LEGACY_OBS_KEYS[role]
        value = deep_get(config, legacy)
        if role == "batch" and value is not None:
            candidates = [str(key) for key in value]
            present = [key for key in candidates if obs is None or key in obs]
            value = present[0] if present else None
        if value is not None:
            values[role], sources[role] = str(value), legacy
        else:
            values[role], sources[role] = default, "default" if default is not None else "unset"
    if obs is not None:
        missing = [f"{role}={values[role]!r} ({sources[role]})" for role in OBS_KEY_ROLES
                   if (role in declared or role in required) and (values[role] is None or values[role] not in obs)]
        if missing:
            columns = list(obs.columns)
            raise KeyError(f"obs keys missing from obs: {', '.join(missing)}; set run.obs_keys to existing columns "
                           f"(obs has {columns[:20]}{' ...' if len(columns) > 20 else ''})")
    return ObsKeys(sources=sources, **values)
