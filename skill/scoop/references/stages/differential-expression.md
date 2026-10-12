# Pseudobulk DE

For patient-level abundance, scSurvival, proportion/UCell, survival or response,
read [outcome association](outcome-association.md) instead of treating the question
as gene-level condition DE. The two worked cases include actual output tables.


Inputs:

- Annotated H5AD.
- Raw counts.
- donor/sample/condition/cell-type metadata.

Scripts:

```bash
scoop de prepare --config runs/<run_id>/config/run.yaml
fastde deseq2 runs/<run_id>/differential-expression/pseudobulk/<cell_type> condition ctrl test
Rscript scripts/r/run_pseudobulk_edger.R <contrast_dir>
```

Default tool:

- FastDE pseudobulk DESeq2-like Wald test. edgeR remains a reference validation path.

Outputs:

- `pseudobulk/counts.tsv`.
- `pseudobulk/metadata.tsv`.
- `contrasts/<contrast>/de_fastde_deseq2.tsv`.

Review triggers:

- Fewer than three biological replicates per group.
- Batch and condition are confounded.
- Rare cell type has too few cells per sample.
