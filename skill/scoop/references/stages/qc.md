# QC by assay

The ordinary scRNA workflow below retains per-sample QC and Scrublet. Imaging
spatial and spot assays use the separate ranges/checklist linked below.


Inputs:

- `artifacts/validated_adata.h5ad` or configured input H5AD.
- `layers["counts"]`.
- `obs["sample_id"]`.

Script:

```bash
scoop qc --config runs/<run_id>/config/run.yaml
```

H5AD writes:

- QC metrics in `obs`.
- `doublet_score`, `doublet_call_scrublet`, `doublet_call`.
- `final_use` updated only for high-confidence failures and doublets.

Sidecar outputs:

- `tables/qc_thresholds_by_sample.tsv`.
- `tables/scrublet_summary_by_sample.tsv`.
- `tables/scrublet_scores.parquet`.

Review triggers:

- Any sample has `fail_fraction > 0.30`.
- QC removes most cells from one condition.
- Ribosomal-high cells are large in number and drive early PCs.

Assay routing: ordinary RNA defaults to `run.assay: scrna`. For Xenium/CosMx and
spot ST, use the separate [spatial QC ranges and checklist](../docs/assay-preprocessing.md).
This route records controls/geometry/coordinates, skips droplet Scrublet, and
never treats unassayed mitochondrial genes as evidence of healthy cells.

Independent libraries use bounded process workers (`qc.scrublet.workers: auto`,
`threads_per_worker: 1`) and count/ID/parameter/version-bound checkpoints under
`qc/scrublet_checkpoints`. The summary includes timing, worker count and resume
state. Use the capture-library key; donors and lineages are not capture libraries.
For exact retained subsets, `--inherit-from parent_qc.h5ad` verifies counts and
axes before inheriting original-library decisions; see [data boundaries](../docs/data-boundaries.md).
