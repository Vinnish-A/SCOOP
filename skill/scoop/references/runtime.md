# Runtime

## Home

The installed product root is `${SCOOP_HOME:-$HOME/.scoop}`. `SCOOP_HOME`, when
set, must be an absolute existing directory. Project data and run artifacts do
not belong in this directory.

Prefer `bin/scoop` and its manifest environment. Unusable CLI support permits a
real API/project adapter; isolate conflicting dependencies outside the shared
release. Record executable/packages, method/reason, inputs/parameters, outputs
and checks; do not claim a managed stage or silently change algorithms.

## Unified environment

For runtime development/release, keep dependencies in `pyproject.toml`'s `runtime`
extra and generate `uv.lock`; no hand-edited lock or secondary environment files.

The unified environment provides core, atlas, OmicVerse, spatial and CCC
capabilities together. Development tools remain a dependency group in the same
`pyproject.toml`, but are excluded from installed releases.

An environment identity includes the frozen dependency lock, Python
implementation/version, platform and machine architecture. Never reuse an
environment whose completion manifest does not match all identity fields.

## Doctor

`scoop doctor --json` reports installed packages, executable stage adapters and each capability separately.
Its `check_kind` and `backend_execution_tested=false` distinguish discovery from
actual execution.
The installer additionally imports every required runtime module and checks the
OmicVerse public APIs before writing the environment completion manifest. Any
managed runtime capability or import-smoke failure rejects the candidate. GPU
availability is hardware-dependent. Before costly analysis, test the actual
backend on representative input: semantics, outputs, memory and timing. Import
success is insufficient; stop hanging probes and use a recorded alternative.

## Installation safety

Build a candidate release in a temporary directory, verify its allowlisted
files and hashes, create/reuse its locked environment, run candidate doctor,
run the full import smoke, then atomically replace the active `current` link.
Refuse to overwrite a non-symlink Agent Skill path. Retain the previous release
for rollback.

For runtime changes, require all of the following before release: regenerate
and check `uv.lock`, perform a cold install in a new temporary SCOOP home, run
the import smoke and CLI help checks, validate release-manifest hashes, and run
the repository test suite. Candidate and production builds from the same tree
must produce the same content ID. Release dependencies stay in `pyproject.toml`;
project-isolated adaptations above do not change the managed release definition.

Release format v3 separates application/Skill content from dependencies. Each
immutable release owns its extracted application and console launchers. A
dependency environment contains only the exported locked dependencies and is
identified by that dependency export plus the resolved Python implementation,
patch version, ABI, platform and architecture. Code/Skill changes reuse the same
dependency environment; they cannot mutate code needed for rollback. Legacy v2
release launchers remain usable as rollback targets. The unified runtime also
includes NetworKit, pinned RCTD-py and official cNMF; RCTD's lookup matrices are cached data,
not another Python environment.

For project Python/API scripts, use `$SCOOP_HOME/bin/python` (or the pinned
`$SCOOP_HOME/current/bin/python`). This launcher uses the same dependency Python
and prepends that release's application modules; it creates no extra environment.
`environment.path` identifies the dependency interpreter, not an installed SCOOP
wheel. Existing legacy environment paths are retained.
