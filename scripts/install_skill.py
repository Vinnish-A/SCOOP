#!/usr/bin/env python3
"""Select an English or original-language SCOOP Skill without rebuilding the runtime."""

import argparse
from datetime import datetime
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
    parser.add_argument("--agent", choices=("codex", "claude"), default="codex")
    parser.add_argument("--agent-home", "--codex-home", dest="agent_home", type=Path,
                        help="Override the selected agent home (--codex-home is a legacy alias)")
    args = parser.parse_args()
    if args.agent == "codex":
        default_agent_home = Path(os.environ.get("CODEX_HOME", Path.home() / ".codex"))
    else:
        default_agent_home = Path(os.environ.get("CLAUDE_CONFIG_DIR", Path.home() / ".claude"))
    agent_home = (args.agent_home or default_agent_home).expanduser().resolve()
    source = ROOT / "skill" / "en" / "scoop" if args.language == "en" else ROOT / "skill" / "scoop"
    if not (source / "SKILL.md").is_file():
        parser.error(f"Missing skill: {source}")
    skill_home = args.home.expanduser().resolve() / "skills" / "scoop"
    revisions = skill_home / "revisions"
    current = skill_home / "current"
    agent_link = agent_home / "skills" / "scoop"
    if args.agent == "codex" and current.exists() and not current.is_symlink():
        parser.error(f"Refusing to overwrite a non-symlink Skill path: {current}")
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
    if agent_link.exists() or agent_link.is_symlink():
        backup = agent_home / "bak" / ("scoop-" + datetime.now().strftime("%Y%m%d-%H%M%S-%f"))
        backup.mkdir(parents=True)
        if agent_link.is_symlink():
            (backup / "scoop").symlink_to(agent_link.resolve())
        else:
            agent_link.rename(backup / "scoop")
        print(f"Previous agent Skill: {backup / 'scoop'}")
    agent_link.parent.mkdir(parents=True, exist_ok=True)
    if args.agent == "codex":
        if current.is_symlink():
            print(f"Previous managed Skill: {current.readlink()}")
        atomic_link(revision, current)
        atomic_link(current, agent_link)
    else:
        # Keep Claude's language choice independent of Codex's managed link.
        atomic_link(revision, agent_link)
    print(f"Installed SCOOP Skill ({args.agent}, {args.language}): {agent_link} -> {revision}")


if __name__ == "__main__":
    main()
