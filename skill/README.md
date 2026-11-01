# SCOOP Skill language editions

- `scoop/`: the original edition, including its existing mixture of Chinese and English.
- `en/scoop/`: the complete English edition, retaining the original structure, evidence requirements, examples, and technical identifiers.

From the repository root, select the edition to install:

```bash
python3 scripts/install_skill.py --language en
python3 scripts/install_skill.py --language zh
# English for Claude Code:
python3 scripts/install_skill.py --agent claude
```

English and Codex are the defaults. `--agent claude` selects Claude Code's
`${CLAUDE_CONFIG_DIR:-~/.claude}/skills/scoop`; `--agent codex` selects
`${CODEX_HOME:-~/.codex}/skills/scoop`. Each agent gets only its selected language.
The selector copies the chosen tree into an immutable content-addressed revision
under `${SCOOP_HOME:-~/.scoop}/skills/scoop/revisions/`. Codex uses the managed
`current` link; Claude is pinned directly to its selected revision, so switching
Claude's language does not change Codex's selection. Both reuse the SCOOP runtime.
Previous agent Skill entries are retained in that agent home's `bak/` directory.
Use `--home PATH` to select a SCOOP home or `--agent-home PATH` to override the
agent home (`--codex-home PATH` remains a compatible alias).
Start a new session and invoke `/scoop` in Claude Code. Its filesystem Skill
format is documented at https://code.claude.com/docs/en/skills.

This is a Skill-only operation. The SCOOP application, dependency environment,
and analysis results are unchanged. The full application installer also defaults
to the English Skill. Use
`python3 scripts/install_scoop.py --skill-language zh` to install the original
edition with the application.

For a fresh Claude-only setup, install the runtime without registering a Codex
Skill, then register the English Skill for Claude Code:

```bash
python3 scripts/install_scoop.py --no-skill-link
python3 scripts/install_skill.py --agent claude
```
