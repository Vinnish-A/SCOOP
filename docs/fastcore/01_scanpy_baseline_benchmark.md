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
| FastCore entry | `fastcore` | `scanpy_legacy` | true | 215.17 s | 1344 MB |

Observed speed ratio:

```text
scanpy_baseline / fastcore_entry = 1.03x
```

This is not a meaningful acceleration. In this benchmark, FastCore only adds
capability planning and manifest/audit output, then runs the same
`scanpy_legacy` backend because OmicVerse adapters are disabled and OmicVerse is
not installed in the Fast environment.

Output consistency:

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

Conclusion:

The current FastCore PR proves the architecture and audit path, not a real
compute acceleration. A fair speedup benchmark requires enabling a validated
OmicVerse/RAPIDS/Rust backend and comparing it against `scanpy_legacy` with the
same preprocessing, HVG, PCA, Harmony, graph, UMAP, and Leiden parameters.
