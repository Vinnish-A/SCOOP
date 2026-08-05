# Runtime

## Home

The installed product root is `${SCOOP_HOME:-$HOME/.scoop}`. `SCOOP_HOME`, when
set, must be an absolute existing directory. Project data and run artifacts do
not belong in this directory.

The managed launcher is `bin/scoop`. It selects an environment using the
requested capability and the active release manifest. Never infer readiness
from a host import or a command found elsewhere on `PATH`.

## Profiles

| Profile | Capabilities |
| --- | --- |
| `core` | validation, QC, FastCore, FastCNMF, FastCNVpy, FastDE, annotation, delivery |
| `atlas` | core plus scvi-tools, scArches/scPoli and atlas reports |
| `spatial` | core plus spatial dependencies |
| `ccc` | core plus FastCCC, LIANA and CellPhoneDB adapters |
| `dev` | tests, lint, type checking and build validation |

An environment identity includes the frozen dependency lock, selected profile,
Python implementation/version, platform and machine architecture. Never reuse
an environment whose completion manifest does not match all identity fields.

## Doctor

`scoop doctor --json` reports each capability separately. A useful check tests
the compatible version and required API, not just module or binary presence.
Optional capability failure must not hide core readiness.

## Installation safety

Build a candidate release in a temporary directory, verify its allowlisted
files and hashes, create/reuse its locked environment, run candidate doctor,
then atomically replace the active `current` link. Refuse to overwrite a
non-symlink Agent Skill path. Retain the previous release for rollback.
