# SCOOP development rules

## Product boundary

- SCOOP is one product, one Skill, and one public `scoop` CLI. Do not rename it
  into an Omni product or introduce a suite manager.
- Install product state under `${SCOOP_HOME:-$HOME/.scoop}`. Keep project data
  and analysis outputs outside that directory.
- Keep orchestration in `skill/scoop/SKILL.md`, detailed biological and
  technical knowledge in `skill/scoop/references/`, deterministic behavior in
  `src/`, and templates or schemas in `skill/scoop/assets/`.

## Compatibility

- Preserve source `docs/` and `markers/` content during migration. Reference
  copies may be reorganized semantically, but no content may be silently
  dropped or summarized.
- Use semantic stage, module, script, output, and reference names. Do not add
  order-prefixed names. Schema versions and biological hierarchy names such as
  `cell_type_lvl1` are allowed.
- Do not add repository-relative `sys.path` changes or assume the current
  working directory. Resolve paths from configuration or explicit arguments.
- Change Python dependencies through `pyproject.toml` and commit `uv.lock`.

## Workflow safety

- Preserve raw counts in `layers["counts"]`; do not use integrated embeddings
  for expression-scale differential testing.
- Every computational stage must create an atomic receipt and expose its review
  state. Never silently change algorithms or dependency profiles.
- scVI integration is followed by clustering, annotation, confidence scoring,
  and the user report. scPoli must verify a current hash-bound approval and may
  supervise only high-confidence labeled cells.
- Agents export and review annotation evidence. Deterministic code validates and
  commits labels to H5AD.

## Release gate

- Run Ruff, type checks, unit tests, CLI contract checks, Skill validation, and
  a synthetic end-to-end test before switching the active release.
- Release only allowlisted Skill files, the built wheel, lock/export files, and
  a content manifest. Reject caches, credentials, data, and temporary files.
- Build candidates in staging, validate there, then atomically switch symlinks.
  Refuse to replace a non-symlink Codex Skill registration and retain the prior
  release for rollback.
