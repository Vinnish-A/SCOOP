from __future__ import annotations

from pathlib import Path

import yaml

from .skill_models import AnnotationSkill


def load_skill(path: str | Path) -> AnnotationSkill:
    path = Path(path)
    skill_yaml = path / "skill.yaml" if path.is_dir() else path
    with skill_yaml.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    if not isinstance(data, dict):
        raise ValueError(f"skill YAML must be a mapping: {skill_yaml}")
    return AnnotationSkill.from_dict(data)


def list_skills(root: str | Path) -> list[AnnotationSkill]:
    root = Path(root)
    if not root.exists():
        return []
    skills = []
    for skill_yaml in sorted(root.glob("*/skill.yaml")):
        skills.append(load_skill(skill_yaml))
    return skills
