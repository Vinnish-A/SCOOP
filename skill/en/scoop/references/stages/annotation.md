# Annotation evidence


Inputs:

- Credible clustering and conventional DEG with tissue marker evidence for broad review.
- For formal fine-subtype review, the accepted biological parent and all evidence
  required by the [annotation workflow](../markers/annotation-workflow.md),
  including batch-handled parent NMF and within-parent clustering/DEG.

Script:

```bash
scoop annotate evidence --config runs/<run_id>/config/run.yaml
```

OmicVerse fast marker option:

```bash
scoop annotate evidence --config runs/<run_id>/config/run.yaml --use-omicverse-cosg
```

The OmicVerse wrapper only exports marker evidence. Delegates may assemble evidence
and propose labels; the main agent personally reads the actual evidence and makes
each final biological decision, which must pass validation before commit. Broad
review does not require NMF; an incomplete fine review remains provisional.

H5AD writes:

- No final labels without the main agent's evidence-based adjudication and validation.

Sidecar outputs:

- `tables/cluster_markers.parquet` (registered with its `label_key`).
- `tables/annotation_evidence_template.tsv`.
- With a COSG option, `tables/cluster_specificity_scores.parquet`; under a
  measurement it scores only genes measured in every sample and lists the rest
  in `tables/cluster_specificity_excluded_genes.tsv`.

Marker DEG reads the registered feature measurement, or `--feature-measurement
<npz>` (which also registers it). Genes measured in every sample go through
scanpy; partially measured genes get an exact tie-corrected Mann–Whitney on
measuring cells only, and need `annotation.default_marker_method: wilcoxon`.
Fewer than 3 measuring cells on either side leaves no row (unassessed). Wilcoxon
output adds `U`, `n_measured_group`, `n_measured_rest`, `log10_pval` and
`log10_pval_adj`, which stay finite where `pvals` underflows.

`annotate review` and `annotate commit` take `--cluster-table` (`obs_name` plus
the cluster-key column) when the input lacks that obs column. Commit also takes
`--cell-layer` and `--output-format obs-parquet`; see the
[annotation workflow](../markers/annotation-workflow.md#structured-decision-contract).

Review triggers:

- Reference and marker evidence disagree.
- New cell type claim.
- Low-confidence annotation used in downstream conclusions.
