# 06. FastCore / OmicVerse-backed 02_core Architecture

FastCore is the new compute engine for the existing `02_core` SOP module. It is
not a new eighth SOP module. The public entry point remains:

```text
scripts/02_core_analysis.py
  -> scsp_agent_sop.core_runner.run_core_pipeline()
```

The runner performs capability planning before execution, selects one backend,
and runs that backend as a whole pipeline. The only non-FastCore fallback is
`scanpy_legacy`, which wraps the original SCOOP normalization, HVG, PCA,
Harmony 2.0, kNN, UMAP, and Leiden sweep implementation.

## Backends

FastCore primary backends:

- `omicverse_cpu`
- `omicverse_cpu_gpu_mixed`
- `omicverse_rust_oom`

Fallback backend:

- `scanpy_legacy`

There is no per-step fallback chain. If an OmicVerse backend is not available or
has not passed quality gates, the planner selects `scanpy_legacy` before the run
starts.

## Capability Planning

`fastcore.runtime.detect_capabilities()` reports optional dependencies without
making imports fail at module import time:

```json
{
  "omicverse_available": true,
  "torch_available": true,
  "cuda_available": true,
  "anndataoom_available": true,
  "rust_backend_available": true,
  "selected_backend": null,
  "fallback_required": false,
  "reasons": []
}
```

`fastcore.backend_plan.plan_fastcore_backend()` then applies the policy:

```text
backed or very large + AnnDataOOM/Rust available -> omicverse_rust_oom
CUDA + torch available                          -> omicverse_cpu_gpu_mixed
External or vendored OmicVerse CPU available   -> omicverse_cpu
otherwise                                       -> scanpy_legacy
```

The default implementation keeps `core.fastcore.enable_omicverse_cpu_backend`
enabled. On ordinary CPU environments this selects the vendored `omicverse_cpu`
backend first; `scanpy_legacy` remains the only fallback when the backend is
disabled or unavailable. Batch correction defaults to Harmony 2.0 through
`harmonypy>=2.0,<3`; Harmony Py Touch is no longer part of the default core
workflow.

The default Fast environment installs `harmonypy>=2.0,<3`, torch/CUDA, CuPy,
and AnnDataOOM. It does not install the external OmicVerse package.
The `omicverse_cpu` backend is a vendored GPL subset of OmicVerse `pp` CPU core code under
`src/fastcore/vendor/omicverse_gpl/`. External OmicVerse validation uses
`environment_omicverse.yml` or the `omicverse` Python extra in a separate
environment.

## Executable Backend Adapters

`omicverse_cpu` is self-contained through the vendored GPL CPU subset.
`omicverse_cpu_gpu_mixed` and `omicverse_rust_oom` are the only accelerated
FastCore adapters. The earlier pure-GPU preprocess path has been removed from
the active FastCore scope.

The mixed backend calls:

```text
ov.settings.cpu_gpu_mixed_init()
ov.pp.preprocess -> ov.pp.scale(use_implicit_centering=True) -> ov.pp.pca
Harmony 2.0 CPU bridge
ov.pp.neighbors(transformer='pyg') -> ov.pp.umap -> ov.pp.leiden
```

The Rust/OOM backend is path based and starts with:

```text
ov.read(input_h5ad, backend='rust')
```

It runs OOM-compatible preprocessing with `shiftlog|pearson`, materializes only
the final core AnnData after PCA, then runs the Harmony/graph/UMAP/Leiden tail
and writes the output H5AD from the runner. The script entry point detects this
backend before reading the input file so the full matrix is not loaded by
`anndata.read_h5ad()` first.

## Stable Output Schema

FastCore must map backend-native results to SCOOP stable keys:

| Stable key | Meaning |
| --- | --- |
| `layers['log1p_norm']` | normalized log expression |
| `obsm['X_pca_biology']` | biology PCA |
| `obsm['X_pca_identity_prebatch']` | identity PCA before batch correction |
| `obsm['X_pca_harmony_identity']` | identity PCA after sample-level correction when available |
| `obsp['connectivities_identity']` | identity kNN graph |
| `obsp['distances_identity']` | identity kNN distances |
| `obsm['X_umap_identity']` | identity UMAP |
| `obsm['X_umap_biology']` | biology UMAP |
| `obs['cluster_identity']` | final selected Leiden cluster |

Large sweep and diagnostic tables stay outside H5AD and are registered in
`adata.uns['file_registry']`.

## Artifacts

Every FastCore run writes:

```text
runs/<run_id>/02_core/fastcore/
  fastcore_manifest.json
  core_quality.json
```

Quality benchmark runs additionally write:

```text
pca_quality.tsv
graph_quality.tsv
cluster_stability.parquet
```

The manifest records backend, fallback status, cell/gene dimensions, HVG/PCA
settings, timing, quality, and external artifact paths.

## Quality Gate

OmicVerse-backed results are accepted only after comparison to a stable
`scanpy_legacy` reference on small and medium datasets:

- PCA explained variance delta
- PC score sign-invariant correlation
- PC subspace cosine
- kNN overlap / Jaccard
- graph density and connected components
- UMAP trustworthiness / local preservation
- Leiden ARI/NMI and seed stability
- wall time, peak RSS, and GPU memory where applicable

Default thresholds are configured under `core.fastcore.quality_gate`.

## Rust / OOM Boundary

The Rust/OOM backend is path-based:

```python
run_omicverse_rust_oom_core(input_h5ad, output_h5ad, cfg, run_root)
```

It must start from `ov.read(path, backend="rust")`, because loading a full
AnnData into memory before calling the backend defeats the out-of-memory design.
`scripts/02_core_analysis.py` therefore performs pre-run planning before any
H5AD read and passes `adata=None` into the runner when `omicverse_rust_oom` is
selected.

## Risks

- OmicVerse is GPL-3.0 licensed. The vendored FastCore backend carries GPL
  provenance and the project declares GPL-3.0-or-later compatibility.
- OmicVerse backend output keys can differ from SCOOP stable keys; FastCore owns
  key mapping.
- GPU torch/CuPy and AnnDataOOM dependencies are optional and must be detected
  before execution.
- Agent/MCP layers must pass paths and manifests, not large matrices.
