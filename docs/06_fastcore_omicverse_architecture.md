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
- `omicverse_gpu_rapids`
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
  "rapids_available": false,
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
CUDA + RAPIDS available                         -> omicverse_gpu_rapids
CUDA + torch available                          -> omicverse_cpu_gpu_mixed
OmicVerse available                             -> omicverse_cpu
otherwise                                       -> scanpy_legacy
```

The initial implementation keeps `core.fastcore.enable_experimental_omicverse_adapters`
disabled by default. This makes the planner deterministic on ordinary CPU
environments and prevents an unvalidated OmicVerse adapter from silently changing
core analysis results. Enabling this flag is a deliberate benchmark/validation
step.

The default Fast environment installs `harmonypy>=2.0,<3` and does not install
OmicVerse. OmicVerse validation uses `environment_omicverse.yml` or the
`omicverse` Python extra in a separate environment.

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

## Risks

- OmicVerse is GPL-3.0 licensed, so distribution and reuse must remain
  license-compliant.
- OmicVerse backend output keys can differ from SCOOP stable keys; FastCore owns
  key mapping.
- GPU, RAPIDS, torch CUDA, and AnnDataOOM dependencies are optional and must be
  detected before execution.
- Agent/MCP layers must pass paths and manifests, not large matrices.
