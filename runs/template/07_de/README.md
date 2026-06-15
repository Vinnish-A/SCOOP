# Pseudobulk DE


Inputs:

- Annotated H5AD.
- Raw counts.
- donor/sample/condition/cell-type metadata.

Scripts:

```bash
python scripts/08_prepare_pseudobulk.py --config runs/<run_id>/config/run.yaml
Rscript r/run_pseudobulk_edger.R <contrast_dir>
```

Default tool:

- edgeR quasi-likelihood model on pseudobulk counts.

Outputs:

- `pseudobulk/counts.tsv`.
- `pseudobulk/metadata.tsv`.
- `contrasts/<contrast>/de_edgeR.tsv`.

Review triggers:

- Fewer than three biological replicates per group.
- Batch and condition are confounded.
- Rare cell type has too few cells per sample.
