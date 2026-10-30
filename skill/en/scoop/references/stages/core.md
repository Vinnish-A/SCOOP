# Core embedding and clustering

Follow the [metadata and batch contract](../docs/design-principles.md#metadata-and-batch-contract):
use the mapped original dataset as the main integration batch key, with verified
technical sub-batches where needed. Bind the actual backend covariates and its
PCA/latent input/output in the receipt; a `study_id` field or mixed UMAP alone does
not establish original-dataset correction or corrected gene expression.

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

Numeric cluster IDs are unsupervised partitions, not biological identities. For
annotation, review the actual clustering and UMAP distances/neighbourhoods with
DEG and tissue markers; two-dimensional UMAP distance alone cannot establish
lineage. A reviewed biological class may combine or split numerical clusters.

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

Validate transfer of a weighted adjacency into the actual clustering backend:
check a nonuniform-edge probe and full-graph total weights/weighted degrees.
A matching input CSR hash does not prove the backend consumed its weights.
NetworKit11.2.2's SciPy-COO overload was observed to insert unit weights even with
`weighted=True`; the corrected explicit `(weights, (rows, columns))` interface
uses float64/uint64 arrays. Record the graph-construction implementation/policy
in checkpoint identity so old unweighted results cannot silently satisfy a fixed
weighted run. Project adapters may use the verified construction when the shared
installed release has not been updated. Expression, PCA/Harmony and the original
weighted CSR remain reusable; only affected clustering descendants change.

A fixed two/four/eight-iteration NetworKit result is an explicit compute budget,
not proof of convergence when its API has no actual termination trace. Preserve
the fast route: use a small coarse grid, independent seeds and budget/stability
or quality comparisons, then deepen only unresolved ranges. An iterative
igraph/leidenalg sensitivity can record real pass gains when needed; do not make
an expensive exhaustive sweep the default just to report convergence everywhere.

Low seed agreement under a fixed budget alone does not establish biological
continuity or prove the inherited representation lost fine-type information.
Inspect actual marker retention and technical structure, then use a bounded
deeper optimization or within-compartment representation comparison when it
resolves a concrete ambiguity. Keep the original valid stages and document which
descendants change; a high neighbour-recall score tests graph search, not biology.

The in-memory core selects observations using validated `final_use` flags before
preprocessing. The Rust/OOM path
applies the same mask through the backed object's lazy row subsetting before
preprocessing; excluded observations never enter fitted representations.
