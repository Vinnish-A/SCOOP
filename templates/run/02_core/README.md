# Core embedding and clustering


Inputs:

- `artifacts/adata_qc.h5ad`.
- `layers["counts"]` and QC flags.

Script:

```bash
python scripts/02_core_analysis.py --config runs/<run_id>/config/run.yaml
```

Default tools:

- Scanpy normalize_total + log1p.
- Scanpy Seurat v3 HVG.
- PCA.
- Torch Harmony on identity PCA.
- Scanpy kNN, UMAP and Leiden.

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
