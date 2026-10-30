# SCOOP Skill language editions

- `scoop/`: the original edition, including its existing mixture of Chinese and English.
- `en/scoop/`: the complete English edition, retaining the original structure, evidence requirements, examples, and technical identifiers.

From the repository root, select the edition to install:

```bash
python3 scripts/install_skill.py --language en
python3 scripts/install_skill.py --language zh
```

The selector copies the chosen tree into an immutable content-addressed revision
under `${SCOOP_HOME:-~/.scoop}/skills/scoop/revisions/`, switches the Skill's
`current` link, and registers it in `${CODEX_HOME:-~/.codex}/skills/scoop`.
It prints the previous revision for rollback and retains all earlier revisions.
Use `--home PATH` or `--codex-home PATH` to override the installation locations.
Start a new session to load the selected edition.

This is a Skill-only operation. The SCOOP application, dependency environment,
and analysis results are unchanged. The full application installer also defaults
to the English Skill. Use
`python3 scripts/install_scoop.py --skill-language zh` to install the original
edition with the application.
