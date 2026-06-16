# Vendored OmicVerse GPL Components

This package contains FastCore's vendored OmicVerse-compatible CPU core.

Upstream:

- Project: OmicVerse
- Version inspected: `omicverse==2.2.3`
- Upstream package path: `omicverse/pp/`
- Relevant files inspected/adapted:
  - `pp/_preprocess.py`
  - `pp/_scale.py`
  - `pp/_pca.py`
  - `pp/_neighbors.py`
  - `pp/_umap.py`
  - `pp/_leiden.py`
- License: GNU General Public License v3.0

FastCore modifications:

- The implementation is reduced to the CPU single-cell core path needed by
  SCOOP `02_core`.
- Optional GPU, Rust/OOM, report, registry, plotting, and lazy pipeline hooks are
  intentionally excluded from this vendored subset.
- Output keys are mapped to SCOOP stable fields.
- Harmony correction uses `harmonypy>=2.0,<3`.
- Leiden uses a single configured resolution instead of the historical SCOOP
  multi-seed sweep; sweep benchmarking remains in `scanpy_legacy`.

Because this vendored backend is derived from GPL-licensed OmicVerse source,
SCOOP/FastCore distribution is GPL-compatible and the project declares
`GPL-3.0-or-later`.
