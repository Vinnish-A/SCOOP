# SCOOP operating guide

SCOOP is the foundation layer of VEngine and is not installed as a separate agent
Skill. This directory holds its operating guide in two editions:

- `scoop/`: the original edition, including its existing mixture of Chinese and English.
- `en/scoop/`: the complete English edition, retaining the original structure, evidence requirements, examples, and technical identifiers.

`scripts/install_scoop.py` ships both editions in every runtime release, as
`${SCOOP_HOME:-~/.scoop}/current/skill/en/` and `.../skill/zh/`. VEngine reads the
edition matching its own language from there, so the guide always describes the
runtime that is actually active, and a rollback restores both together. Install
VEngine from the V-Skills repository; it is the only agent Skill for this system.
