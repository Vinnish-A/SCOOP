# FastCore versus Scanpy Legacy Baseline

Benchmark date: 2026-06-16

Environment:

- Python: `.venv-scoop-fast`, CPython 3.11.12
- Batch correction: `harmonypy==2.0.0`
- OmicVerse: not installed in the Fast environment
- Torch / Harmony-PyTorch: not installed

Input:

- `h5ad/canonical/quick_test/public_O_GSE154795_10samples_1000cells.h5ad`
- Shape: 10,000 cells x 31,402 genes
- Counts layer: `counts`

Configuration:

- HVG: 3,000 genes
- PCA: 50 components
- kNN: 15 neighbors
- UMAP: `min_dist=0.3`
- Leiden sweep: 8 resolutions x 5 seeds
- Batch correction: Harmony 2.0, `max_iter_harmony=20`

Results:

| Run | Engine | Selected backend | Fallback used | Wall time | Peak RSS |
| --- | --- | --- | --- | ---: | ---: |
| scanpy baseline | `scanpy_legacy` | `scanpy_legacy` | false | 222.10 s | 1399 MB |
| FastCore entry before vendoring | `fastcore` | `scanpy_legacy` | true | 215.17 s | 1344 MB |
| FastCore vendored OmicVerse CPU | `fastcore` | `omicverse_cpu` | false | 52.61 s | 1469 MB |

Observed speed ratio:

```text
scanpy_baseline / fastcore_entry_before_vendoring = 1.03x
scanpy_baseline / fastcore_vendored_omicverse_cpu = 4.22x
```

The first FastCore entry was not a meaningful acceleration because it only added
capability planning and manifest/audit output, then ran the same
`scanpy_legacy` backend. The vendored OmicVerse CPU backend is a real compute
change and removes the historical 8-resolution x 5-seed Leiden sweep from the
default FastCore path, matching OmicVerse's single-resolution core workflow.

Fallback output consistency before vendoring:

| Metric | Value |
| --- | ---: |
| Cluster ARI | 1.0 |
| Cluster NMI | 1.0 |
| `X_pca_biology` MAE | 0.0 |
| `X_pca_identity_prebatch` MAE | 0.0 |
| `X_pca_harmony_identity` MAE | 0.0 |
| `X_umap_identity` MAE | 0.0 |
| `X_umap_biology` MAE | 0.0 |

Main wall-time steps from the FastCore entry run:

| Step | Seconds |
| --- | ---: |
| `leiden_sweep` | 155.73 |
| `neighbors_umap_identity` | 40.86 |
| `neighbors_umap_biology` | 9.96 |
| `pca_biology` | 1.21 |
| `pca_identity_prebatch` | 1.14 |
| `select_hvg_biology` | 1.00 |
| `harmony2` | 0.99 |

Vendored OmicVerse CPU 10k wall-time steps:

| Step | Seconds |
| --- | ---: |
| `neighbors_umap_single` | 40.54 |
| `leiden_single` | 4.69 |
| `shiftlog_normalize` | 1.37 |
| `pca_covariance_eigh` | 1.20 |
| `hvg_seurat` | 0.99 |
| `harmony2` | 0.99 |
| `scale_hvg` | 0.48 |

OmicVerse reference similarity on a 2k-cell subset:

| Metric | Value |
| --- | ---: |
| HVG Jaccard | 1.0 |
| PCA 50-PC subspace cosine | 1.0 |
| Harmony 50-PC subspace cosine | 0.99999999995 |
| kNN graph edge Jaccard | 1.0 |
| Cluster ARI | 0.976 |
| Cluster NMI | 0.964 |

External backend smoke tests:

| Backend | Environment | Shape | Wall time | Peak RSS | Status |
| --- | --- | ---: | ---: | ---: | --- |
| `omicverse_cpu_gpu_mixed` | `.venv-scoop-omicverse` | 120 x 300 | 13.91 s | 2075 MB | passed |
| `omicverse_rust_oom` | `.venv-scoop-omicverse` + `anndataoom==0.1.8` | 120 x 300 | 20.37 s | 1199 MB | passed |
| `omicverse_gpu_rapids` | not installed in current env | n/a | n/a | n/a | adapter covered by mock test; real RAPIDS env still required |

The mixed smoke confirmed the OmicVerse torch/pyg PCA-neighbors-UMAP path and
the FastCore Harmony 2.0 CPU bridge. The Rust/OOM smoke confirmed that
`scripts/02_core_analysis.py` avoids pre-reading the H5AD and starts through
`ov.read(..., backend='rust')`; graph/UMAP/Leiden still run after the final
minimal AnnData materialization.

Conclusion:

The vendored OmicVerse CPU backend gives a measured `4.22x` speedup over the
historical fallback on the 10k public quick test, while matching the external
OmicVerse CPU reference closely on the tested 2k subset. The remaining major
bottleneck is still UMAP/neighbors.
