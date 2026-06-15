# QC, Scrublet and ambient flags


Inputs:

- `artifacts/validated_adata.h5ad` or configured input H5AD.
- `layers["counts"]`.
- `obs["sample_id"]`.

Script:

```bash
python scripts/01_qc_scrublet.py --config runs/<run_id>/config/run.yaml
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
