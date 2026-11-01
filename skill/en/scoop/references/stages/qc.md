# QC by assay

The ordinary scRNA workflow below retains per-sample QC and Scrublet. Imaging
spatial and spot assays use the separate ranges/checklist linked below.


Inputs:

- `artifacts/validated_adata.h5ad` or configured input H5AD.
- `layers["counts"]`.
- `obs["sample_id"]`, or the declared `run.obs_keys.sample` (resolved once; a declared key must exist).

Script:

```bash
scoop qc --config runs/<run_id>/config/run.yaml
```

H5AD writes:

- QC metrics in `obs`. `pct_counts_mt/ribo/hb` is NaN and `mt/ribo/hb_unassessed`
  True when the axis has no flagged gene or, with a registered
  [feature measurement](../docs/data-boundaries.md#observability-positive-evidence-and-input-versions),
  the cell's sample measured none. Unassessed is not 0%; NaN never crosses a threshold.
  MT/ribo/HB gene flags read symbols from `var[qc.symbol_key]` when set; a
  `qc.symbol_key` that is not a var column is an error, not a fallback to var_names.
- `doublet_score`, `doublet_call_scrublet`, `doublet_call`.
- `final_use` updated only for high-confidence failures and doublets.
- QC-owned obs columns already in the input (`qc_*`, `doublet_*`, `pct_counts_*`,
  `log10_*`, `total_counts`, `n_genes_by_counts`, `*_unassessed`, `final_use`, ...)
  are renamed to `prior_<name>` before QC recomputes them (an older `prior_<name>`
  is replaced) and listed in the decision log (S026). The new `final_use` is ANDed
  with `prior_final_use`: QC never re-admits a cell an earlier step excluded.
- `--output-format obs-parquet --output <path>.parquet` writes only the QC obs
  columns (and the `prior_*` they displaced) plus a JSON receipt instead of the H5AD.

Sidecar outputs:

- `qc/tables/qc_thresholds_by_sample.tsv`.
- `qc/tables/scrublet_summary_by_sample.tsv`.
- `qc/tables/scrublet_scores.parquet`.

Review triggers:

- Any sample has `fail_fraction > 0.30`.
- QC removes most cells from one condition.
- Ribosomal-high cells are large in number and drive early PCs.

Assay routing: ordinary RNA defaults to `run.assay: scrna`. For Xenium/CosMx and
spot ST, use the separate [spatial QC ranges and checklist](../docs/assay-preprocessing.md).
This route records controls/geometry/coordinates, skips droplet Scrublet, and
never treats unassayed mitochondrial genes as evidence of healthy cells.
Spatial QC, spatial core, RCTD and `spatial build-reference` take gene features from
`var.is_gene_expression`, else `feature_types` matched against `qc.gene_feature_types`
(default `[gene expression, gene]`).

Independent libraries use bounded process workers (`qc.scrublet.workers: auto`,
`threads_per_worker: 1`) and count/ID/parameter/version-bound checkpoints under
`qc/scrublet_checkpoints`. The summary includes timing, worker count and resume
state. Use the capture-library key; donors and lineages are not capture libraries.
For exact retained subsets, `--inherit-from parent_qc.h5ad` verifies counts and
axes before inheriting original-library decisions; see [data boundaries](../docs/data-boundaries.md).


## Reused QC and lineage-aware reference admission

Resolve effective QC from actual executed capture/library receipts and exact-UID outputs. A stale summary such as `not_run` is not evidence that an available completed assessment never ran; retain that provenance discrepancy rather than rerunning valid QC or recoding unknown calls as singlets.

A universal detected-gene cutoff can remove an entire genuinely low-RNA lineage, such as erythroid RNA. Review native count depth, available-axis gene complexity, source singlet/QC evidence and independent multi-gene lineage support together. Any lineage-specific exception is an explicit analysis/reference admission rule with denominators and uncertainty; it does not establish intact cells, normality or an automatically changed whole-model mask.
