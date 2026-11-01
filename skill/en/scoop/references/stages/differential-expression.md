# Gene-level differential expression

For patient-level abundance, scSurvival, proportion/UCell, survival or response,
read [outcome association](outcome-association.md) instead of treating the question
as gene-level condition DE. The two worked cases include actual output tables.


## Cluster annotation markers

Use the [metadata and batch contract](../docs/design-principles.md#metadata-and-batch-contract).
After verifying the clustering, input expression and measurement scope, perform
ordinary gene-level cluster, subcluster and merged final-label comparisons on the
actual expression layer. Batch-stratified tests are optional diagnostics for a
specific source or measurement concern, not a prerequisite for annotation.
Include every source, including collection-derived datasets, when its assay and
gene observability support the contrast. Record actual denominators, effects,
source composition and the tested-gene/BH universe; missing contrasts and
unmeasured genes remain unassessed. Preserve the semantics of counts and
normalized-only inputs; do not invent counts. Interpret source confounding and
measurement differences explicitly. These cell-level contrasts are descriptive
annotation evidence, not independent patient or condition inference.
With a declared feature measurement, `annotate evidence` tests partially measured
genes on measuring cells only and exports `n_measured_*`, `U` and `log10_pval(_adj)`;
see [annotation evidence](annotation.md).

## Patient/condition DE

Use genuine biological replicates and a declared donor/sample/technical-batch
and condition design. Preserve paired patients when known; an unresolved sample
is not an invented donor, and a fully confounded batch/condition effect is not
identifiable. Cell-level marker P values do not provide patient-level replication.

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

FastDE fits `~ condition` only (no batch or donor terms). It builds the IRLS working
response from the floored mean, as DESeq2's `fitBeta` does, so a gene that is zero in
one group converges to a finite log2FC instead of growing with `--max-iter`. The
`betaConv` column flags genes that missed `--beta-tol` or have |beta| > 30 (log2);
the latter get NA `stat`, `pvalue` and `padj` and are left out of the BH universe.
The manifest records `n_beta_not_converged` and `n_beta_out_of_bounds`.

Outputs:

- `pseudobulk/counts.tsv`.
- `pseudobulk/metadata.tsv`.
- `contrasts/<contrast>/de_fastde_deseq2.tsv`.

Review triggers:

- Fewer than three biological replicates per group.
- Batch and condition are confounded.
- Rare cell type has too few cells per sample.

`de prepare` resolves `run.obs_keys` (sample required). When `annotation_status`
exists it aggregates only `accepted`/`labeled` cells; `review_required` (including
reviewed unresolved), `excluded_technical` and `unassigned` cells never enter a pseudobulk. The metadata carries
the donor, batch, study and `keys.batch_candidates` columns present, for the design.
