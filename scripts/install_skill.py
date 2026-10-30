#!/usr/bin/env python3
"""Select an English or original-language SCOOP Skill without rebuilding the runtime."""

import argparse
import hashlib
import os
from pathlib import Path
import shutil
import tempfile
import uuid

ROOT = Path(__file__).resolve().parents[1]


def atomic_link(target, link):
    temporary = link.with_name(link.name + "." + uuid.uuid4().hex + ".tmp")
    temporary.symlink_to(target)
    os.replace(temporary, link)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--language", choices=("en", "zh"), default="en")
    parser.add_argument("--home", type=Path, default=Path(os.environ.get("SCOOP_HOME", Path.home() / ".scoop")))
    parser.add_argument("--codex-home", type=Path, default=Path(os.environ.get("CODEX_HOME", Path.home() / ".codex")))
    args = parser.parse_args()
    source = ROOT / "skill" / "en" / "scoop" if args.language == "en" else ROOT / "skill" / "scoop"
    if not (source / "SKILL.md").is_file():
        parser.error(f"Missing skill: {source}")
    skill_home = args.home.expanduser().resolve() / "skills" / "scoop"
    revisions = skill_home / "revisions"
    current = skill_home / "current"
    codex_link = args.codex_home.expanduser().resolve() / "skills" / "scoop"
    for link in (current, codex_link):
        if link.exists() and not link.is_symlink():
            parser.error(f"Refusing to overwrite a non-symlink Skill path: {link}")
    digest = hashlib.sha256()
    for path in sorted(p for p in source.rglob("*") if p.is_file()):
        digest.update(path.relative_to(source).as_posix().encode() + b"\0")
        digest.update(path.read_bytes())
    revision = revisions / (args.language + "_" + digest.hexdigest()[:20])
    revisions.mkdir(parents=True, exist_ok=True)
    if not revision.exists():
        with tempfile.TemporaryDirectory(prefix=".skill-", dir=revisions) as temporary:
            staging = Path(temporary) / "skill"
            shutil.copytree(source, staging)
            staging.rename(revision)
    if current.is_symlink():
        print(f"Previous Skill: {current.readlink()}")
    codex_link.parent.mkdir(parents=True, exist_ok=True)
    atomic_link(revision, current)
    atomic_link(current, codex_link)
    print(f"Installed SCOOP Skill ({args.language}): {revision}")


if __name__ == "__main__":
    main()
