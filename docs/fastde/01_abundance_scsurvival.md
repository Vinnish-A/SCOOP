# FastDE Abundance

`fastde abundance` performs sample-level differential cell abundance and outcome
association. It is not gene-level DE and it never uses per-cell p values.

Input can be either:

- an H5AD with `obs[sample_key]` and `obs[celltype_key]`;
- a sample-by-cell-type count matrix TSV.

The module constructs:

- `Y`: sample x cell-type counts;
- `P`: sample x cell-type proportions;
- `N`: total cells per sample;
- `M`: sample metadata;
- `X`: transformed abundance features plus optional covariates.

## Modes

- `survival`: scSurvival-like Cox survival association from abundance features.
- `binary`: scSurvival-like binary outcome classifier.
- `multiclass`: subtype/class association from abundance features.
- `condition`: alias of binary condition comparison.
- `continuous`: association with a continuous phenotype.

Default abundance transform is CLR:

```text
P = (Y + 0.5) / row_sum(Y + 0.5)
Z = log(P) - row_mean(log(P))
```

## Examples

Survival:

```bash
fastde abundance \
  --mode survival \
  --input-h5ad runs/<run_id>/artifacts/adata_annotation_committed.h5ad \
  --sample-key sample_id \
  --celltype-key cell_type_lvl3 \
  --metadata runs/<run_id>/config/sample_metadata.tsv \
  --time-col OS_time \
  --event-col OS_event \
  --covariates age,sex,batch \
  --output-dir runs/<run_id>/07_de/abundance_survival
```

Binary:

```bash
fastde abundance \
  --mode binary \
  --input-h5ad runs/<run_id>/artifacts/adata_annotation_committed.h5ad \
  --sample-key sample_id \
  --celltype-key cell_type_lvl3 \
  --metadata runs/<run_id>/config/sample_metadata.tsv \
  --label-col responder \
  --positive-label response \
  --negative-label non_response \
  --covariates age,sex,batch \
  --output-dir runs/<run_id>/07_de/abundance_binary
```

Multiclass:

```bash
fastde abundance \
  --mode multiclass \
  --input-h5ad runs/<run_id>/artifacts/adata_annotation_committed.h5ad \
  --sample-key sample_id \
  --celltype-key cell_type_lvl3 \
  --metadata runs/<run_id>/config/sample_metadata.tsv \
  --label-col subtype \
  --reference-level control \
  --covariates batch,sex \
  --output-dir runs/<run_id>/07_de/abundance_multiclass
```

Direct matrix:

```bash
fastde abundance \
  --mode survival \
  --counts sample_by_celltype_counts.tsv \
  --metadata sample_metadata.tsv \
  --time-col OS_time \
  --event-col OS_event \
  --output-dir out
```

## Outputs

Every run writes external artifacts:

- `sample_by_celltype_counts.tsv`
- `sample_by_celltype_proportions.tsv`
- `abundance_metadata_used.tsv`
- `abundance_<mode>_results.tsv`
- `abundance_<mode>_predictions.tsv`
- `abundance_<mode>_metrics.json`
- `abundance_manifest.json`

No long abundance result tables are stored inside H5AD.
