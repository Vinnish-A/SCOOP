# Development contract

## Ownership

- Keep one public `scoop` CLI and one `scoop` Skill.
- Keep orchestration in `SKILL.md`, detailed knowledge in references, reusable
  deterministic behavior in runtime code, and field contracts in schemas.
- Give each rule one owner. Do not duplicate defaults between config and a run
  template; generate the latter from the former.
- Preserve Fast engine packages as compute backends. Do not let them decide
  biological labels.

## Naming

- Use semantic stage and artifact names; do not use ordering prefixes.
- Use `snake_case` for Python modules and `kebab-case` for Skill/reference paths.
- Version/schema numerals and biological hierarchy fields such as
  `cell_type_lvl1` are not ordering prefixes and remain valid.

## Dependencies and paths

- Change dependencies only through `pyproject.toml` and `uv.lock`.
- Do not add `sys.path` mutations or repository-CWD assumptions.
- Resolve project paths relative to the config file or as explicit absolute paths.
- Keep optional profiles explicit; never silently install or switch a profile.

## Contracts and compatibility

- Version every configuration, receipt, H5AD, report and approval schema.
- Add a migration or compatibility test for any schema change.
- Treat a CLI help string, reference command and implementation as one tested contract.
- Keep source provenance outside content identity. Release identity depends only
  on allowlisted file content and declared compatibility inputs.

## Quality gates

- Run Ruff, type checks, unit tests, contract-negative tests and a synthetic
  end-to-end workflow before release.
- Test every declared engine command with `--help` or a capability probe.
- Compare semantic outputs against the frozen pre-migration baseline.
- Keep performance benchmarks separate from correctness tests.
- Require source, version, hash and license for vendored code and external assets.
- Reject untracked, cache, data, credential and temporary files from releases.
- Forward-test standard, tumor and atlas requests in isolated contexts.
