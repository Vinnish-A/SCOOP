# Core embedding and clustering


Inputs:

- `artifacts/adata_qc.h5ad`.
- `layers["counts"]` and QC flags.

Script:

```bash
scoop core --config runs/<run_id>/config/run.yaml
```

Default tools:

- Scanpy normalize_total + log1p (`log1p_norm` at `core.normalize_total_target_sum`, every backend).
- Scanpy Seurat v3 HVG.
- PCA.
- Harmony on identity PCA: harmonypy 2.x on CPU, for the declared `run.obs_keys.batch`
  (else the present `keys.batch_candidates`); `fastcore/fastcore_manifest.json` records what ran.
- Scanpy kNN, UMAP and Leiden (`core.leiden_search.selection`: `coarsest_stable` or `finest_stable`).

H5AD writes:

- `layers["log1p_norm"]`.
- `obsm["X_pca_biology"]`.
- `obsm["X_pca_harmony_identity"]`.
- `obsm["X_umap_identity"]` and `obsm["X_umap_biology"]`.
- `obs["cluster_identity"]`, and `obsm["cluster_identity_sweep_labels"]` (one labelling per swept resolution).

Sidecar outputs:

- `tables/hvg_biology.parquet`.
- `tables/hvg_identity.parquet`.
- `tables/leiden_sweep.parquet`.
- `tables/cluster_stability.parquet`.

Review triggers:

- Harmony unavailable (`batch_correction.status: fallback_uncorrected`) or over-corrects markers.
- Selected clustering is unstable across seeds.
- Clusters are dominated by ribosomal/stress/cell-cycle markers.
