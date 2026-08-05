# Runtime

## Home

The installed product root is `${SCOOP_HOME:-$HOME/.scoop}`. `SCOOP_HOME`, when
set, must be an absolute existing directory. Project data and run artifacts do
not belong in this directory.

The managed launcher is `bin/scoop`. It invokes the single unified environment
recorded by the active release manifest. Never infer readiness from a host
import or a command found elsewhere on `PATH`.

## Unified environment

Maintain runtime dependencies only in `pyproject.toml` under the `runtime`
extra. Generate `uv.lock`; never edit it manually. Do not maintain Conda or
secondary environment files.

The unified environment provides core, atlas, OmicVerse, spatial and CCC
capabilities together. Development tools remain a dependency group in the same
`pyproject.toml`, but are excluded from installed releases.

An environment identity includes the frozen dependency lock, Python
implementation/version, platform and machine architecture. Never reuse an
environment whose completion manifest does not match all identity fields.

## Doctor

`scoop doctor --json` reports installed packages and each capability separately.
The installer additionally imports every required runtime module and checks the
OmicVerse public APIs before writing the environment completion manifest. Any
managed runtime capability or import-smoke failure rejects the candidate. GPU
availability remains hardware-dependent and is reported separately.

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
must produce the same content ID. Never add a second environment definition to
work around a dependency conflict; resolve the conflict in `pyproject.toml`.
