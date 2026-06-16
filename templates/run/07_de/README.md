# Pseudobulk DE


Inputs:

- Annotated H5AD.
- Raw counts.
- donor/sample/condition/cell-type metadata.

Scripts:

```bash
python scripts/08_prepare_pseudobulk.py --config runs/<run_id>/config/run.yaml
fastde deseq2 runs/<run_id>/07_de/pseudobulk/<cell_type> condition ctrl test
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
