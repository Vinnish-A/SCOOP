# Core embedding and clustering


Inputs:

- `artifacts/adata_qc.h5ad`.
- `layers["counts"]` and QC flags.

Script:

```bash
scoop core --config runs/<run_id>/config/run.yaml
```

Default tools:

- Scanpy normalize_total + log1p.
- Scanpy Seurat v3 HVG.
- PCA.
- Torch Harmony on identity PCA.
- Scanpy kNN/UMAP and a concurrent NetworKit Leiden resolution/seed sweep.

H5AD writes:

- `layers["log1p_norm"]`.
- `obsm["X_pca_biology"]`.
- `obsm["X_pca_harmony_identity"]`.
- `obsm["X_umap_identity"]` and `obsm["X_umap_biology"]`.
- `obs["cluster_identity"]`.

Sidecar outputs:

- `tables/hvg_biology.parquet`.
- `tables/hvg_identity.parquet`.
- `tables/leiden_sweep.parquet`.
- `tables/cluster_stability.parquet`.

Review triggers:

- Harmony unavailable or over-corrects markers.
- Selected clustering is unstable across seeds.
- Clusters are dominated by ribosomal/stress/cell-cycle markers.

Ordinary scRNA retains its primary method. See [assay-specific preprocessing](../docs/assay-preprocessing.md)
for the recorded small-lineage LOESS recovery and spatial baseline. Spatial
`run.assay` values bypass the RNA FastCore planner; controls never enter PCA.

Managed core backends all use the real shared Leiden resolution/seed sweep.
New defaults run up to four concurrent jobs with one Leiden thread each and a
fixed two-iteration budget. Concurrency is bounded by available CPU/load and
memory; the grid/seeds and configured coarse/refine stability policy are kept.
Override `core.leiden_search.iterations` for a deeper sensitivity run. Per-job
checkpoint metadata records actual workers/threads, runtime version, objective,
iterations and elapsed time. A single seed has unassessed (NaN) stability, never
fabricated ARI/NMI=1. Existing runs/configs are not rewritten.

The in-memory core selects observations using validated `final_use` flags before
preprocessing. The Rust/OOM path
applies the same mask through the backed object's lazy row subsetting before
preprocessing; excluded observations never enter fitted representations.
