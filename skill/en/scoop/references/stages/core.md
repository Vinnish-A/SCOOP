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
- Harmony on the identity PCA: harmonypy 2.x on CPU (`harmonypy.run_harmony`), not Torch Harmony.
- Scanpy kNN/UMAP and a concurrent NetworKit Leiden resolution/seed sweep.

H5AD writes:

- `layers["log1p_norm"]`.
- `obsm["X_pca_biology"]`.
- `obsm["X_pca_harmony_identity"]`.
- `obsm["X_umap_identity"]` and `obsm["X_umap_biology"]`.
- `obs["cluster_identity"]`.

Every RNA backend delivers `log1p_norm` as plain normalize_total at
`core.normalize_total_target_sum` (default 10000) then log1p: `scanpy_legacy` and
`fastcore_cpu` write it themselves, and the core runner writes it after
`fastcore_mixed`/`fastcore_oom` (replacing any input layer of that name), so
`scoop programs` can read it whichever backend the planner picks (S007). The
OmicVerse/`fastcore_cpu` shiftlog matrix (`core.fastcore.omicverse.target_sum`)
stays internal to HVG scaling and PCA. The library size is the counts-layer row
sum unless `core.size_factor_key` names an obs column of per-cell totals (for
example full-gene `total_counts` when the matrix holds a gene subset, S044); a
column smaller than the row sums is refused. `fastcore_manifest.json` records
`log1p_norm` (target sum, size-factor key, writer). Before reporting success the
runner checks the output contract below plus the counts layer and fails naming
anything missing.

The planner sizes an unloaded H5AD from the obs index and the nnz of
`qc.counts_layer` (else X), so a file with counts only in layers plans normally.

Core Harmony batch rule: a declared `run.obs_keys.batch` is the one key Harmony
corrects for. Without it the pre-R10 rule remains: every present
`keys.batch_candidates` column (default the sample key), recorded with source
`keys.batch_candidates` because it is a guess, not a declared batch. Keys with one
level are dropped. `fastcore_manifest.json` `batch_correction` records method,
key source, requested and used keys, dropped keys and the reason when none were
used, plus harmonypy's version, iterations (`n_iter`, `reached_max_iter`) and
objective trace (S046). `status` is `corrected`, `not_run` (no two-level key),
`not_requested` (`method: none`) or `fallback_uncorrected`: harmonypy did not
import, `X_pca_harmony_identity` is the uncorrected PCA, and the stage requests
review. `harmony2_used` is true only for `corrected`.

All backends honour `core.neighbors_n_pcs` (OmicVerse backends also read
`core.fastcore.omicverse.neighbors.n_pcs`); unset, `scanpy_legacy` and
`fastcore_cpu` use every component.

With a registered [feature measurement](../docs/data-boundaries.md#observability-positive-evidence-and-input-versions),
`scanpy_legacy` and `fastcore_cpu` rank each HVG batch only on genes all its
samples measured (adding `n_batches_measured` when any gene is partial).
`scanpy_legacy` programme scores draw set and control genes from genes measured
in every sample (`uns["program_score_genes"]`). Other backends record
`measurement: not_honoured` in `fastcore_manifest.json` and request review.
Rescoring replaces the cell-cycle columns instead of duplicating them.

Numeric cluster IDs are unsupervised partitions, not biological identities. For
annotation, review the actual clustering and UMAP distances/neighbourhoods with
DEG and tissue markers; two-dimensional UMAP distance alone cannot establish
lineage. A reviewed biological class may combine or split numerical clusters.

Sidecar outputs:

- `tables/hvg_biology.parquet`.
- `tables/hvg_identity.parquet`.
- `tables/leiden_sweep.parquet`.
- `tables/cluster_stability.parquet` (with `selection_policy`, `eligible`, `chosen`).

Leiden auto-selection takes, among resolutions whose median seed ARI and cluster
count pass `core.leiden_search.{min_seed_ari,min_clusters,max_clusters}`, the
coarsest (`core.leiden_search.selection: coarsest_stable`, default) or the finest
(`finest_stable`, B07); coarse-to-fine refines around that choice. The policy,
rule, chosen resolution/seed and candidate table are in `fastcore_manifest.json`
`clustering` and `uns["leiden_selection"]["cluster_identity"]`.
`obsm["cluster_identity_sweep_labels"]` keeps one labelling per swept resolution
(seed 0 where it ran) so a reviewed resolution can be chosen explicitly instead.

Review triggers:

- Harmony unavailable (`batch_correction.status: fallback_uncorrected`) or over-corrects markers.
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
