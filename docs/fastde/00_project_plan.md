# FastDE

FastDE separates three analyses that should not be mixed:

- marker genes: "which genes define this cell group?"
- condition DE: "what changed in the same cell group after condition changes?"
- abundance: "which cell-type or cell-state proportions vary with sample-level
  design or outcome?"

## Marker Genes

Default method: COSG-like cosine score.

Input is a cell-by-gene expression matrix and a group label. For every group,
FastDE scores each gene by cosine similarity between the gene expression vector
and the binary group membership vector. The output is a ranked marker table with
`group`, `rank`, `gene`, `score`, `pct_in`, `pct_out`, `mean_in`, `mean_out`,
and `mean_diff`.

Fallback method: Wilcoxon rank-sum. It is slower and intended for confirmation
or compatibility checks, not as the default marker engine.

Sparse acceleration:

- Dense reference materializes the full `cells x genes` matrix and is only safe
  for small subsets.
- FastDE COSG detects sparse input and uses a CSC column scan. It computes
  `group x gene` summaries (`sum`, `nnz`, total sum, and L2 norm) without
  densifying the full expression matrix.
- Peak memory is therefore driven by the input sparse matrix plus small
  `groups x genes` arrays, not by `cells x genes` dense expansion.

CLI:

```bash
fastde markers \
  --input annotated.h5ad \
  --groupby cell_type_lvl3 \
  --layer log1p_norm \
  --output-dir markers \
  --method cosg
```

Benchmark command:

```bash
PYTHONPATH=src python scripts/fastde/benchmark_cosg.py \
  --input .scoop_local/data/h5ad/canonical/quick_test/public_O_GSE154795_24samples_3000cells_balanced.h5ad \
  --groupby sample_id \
  --layer counts \
  --output-dir tmp/fastde_cosg_benchmark_public \
  --full-sparse
```

Current local real-data result:

| Dataset | Shape | Mode | Algorithm | Total | Peak RSS | Dense consistency |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| public quick test subset | 6k x 5k | dense reference | 5.54s | 9.74s | 4.05GB | reference |
| public quick test subset | 6k x 5k | sparse CSC scan | 0.10s | 4.64s | 3.84GB | 1.00 top50/rank |
| public quick test full | 72k x 40,791 | sparse CSC scan | 2.10s | 29.04s | 3.65GB | subset-validated |
| private quick test subset | 6k x 5k | dense reference | 5.12s | 17.87s | 5.50GB | reference |
| private quick test subset | 6k x 5k | sparse CSC scan | 0.11s | 1.99s | 5.36GB | 1.00 top50/rank |
| private quick test full | 72k x 36,581 | sparse CSC scan | 2.65s | 42.43s | 5.02GB | subset-validated |

The remaining full-data total time is dominated by H5AD I/O. The marker scoring
itself is already a few seconds on these 72k-cell objects.

## Condition DE

Default method: Python DESeq2-like pseudobulk NB Wald.

FastDE keeps SCOOP's pseudobulk design. It does not run per-cell condition tests
and does not use integrated expression. The input is one cell-type pseudobulk
directory containing `counts.tsv` and `metadata.tsv`.

Implemented DESeq2 concepts:

- median-ratio size factors.
- design-aware gene-wise negative-binomial dispersion estimates.
- parametric dispersion trend fitting.
- dispersion MAP shrinkage with a log-normal prior.
- dispersion outlier handling.
- log-link GLM with size-factor offset.
- Wald test for the condition coefficient.
- Benjamini-Hochberg adjusted p values.

Not yet implemented:

- LFC shrinkage. This matches DESeq2's default `betaPrior=FALSE` result path,
  but does not replace `lfcShrink`.
- complex multi-factor designs.

CLI:

```bash
fastde deseq2 runs/<run_id>/07_de/pseudobulk/<cell_type> condition ctrl test
```

The output table uses DESeq2-style columns: `gene`, `baseMean`,
`log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`, and `dispersion`.
It also writes dispersion diagnostics: `dispGeneEst`, `dispFit`, `dispMAP`,
and `dispOutlier`.

## DESeq2 Reference Benchmark

Benchmark command:

```bash
PYTHONPATH=src python scripts/fastde/benchmark_deseq2.py \
  --output-dir tmp/fastde_deseq2_benchmark_30k \
  --n-genes 30000 \
  --n-per-group 6 \
  --n-de 600
```

Current local result after installing R `DESeq2`:

| Fixture | FastDE | R DESeq2 | Speedup | log2FC Spearman | Sign agreement | Top100 overlap |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 10k genes, 12 samples | 6.34s | 9.86s | 1.56x | 0.999997 | 0.9996 | 0.98 |
| 30k genes, 12 samples | 18.55s | 19.68s | 1.06x | 0.999994 | 0.9993 | 1.00 |

The high log2FC agreement is expected because both methods use size-factor
offsets and negative-binomial count models. The complete dispersion path is much
slower than the earlier moments-only prototype, but it raises top-hit agreement
from roughly 0.70-0.73 to 0.98-1.00 on these fixtures.

## Abundance

Default method: scSurvival-style multiple-instance learning over sample bags.

`fastde abundance` compares sample-level cell-type or cell-state abundance
across survival, binary, multiclass, continuous, or condition designs. It is not
gene-level DE and does not use per-cell p values. The default backend inherits
the key scSurvival architecture: every sample is a bag of cell instances, a
shared instance encoder maps cells/states into latent features, gated attention
pools each bag into a sample embedding, and the task-specific head predicts
survival risk, binary label, multiclass subtype, or continuous phenotype.

Supported modes:

- `survival`: Cox-style survival association from cell-type proportions.
- `binary`: binary outcome classification and per-cell-type association table.
- `multiclass`: subtype/class association with reference-level contrasts.
- `condition`: binary condition comparison alias.
- `continuous`: continuous phenotype association.
- `--survival-loss cox`: original scSurvival-style Cox partial likelihood.
- `--survival-loss cox_rank`: pairwise ranking loss on comparable survival
  pairs.
- `--survival-loss cox_plus_rank`: Cox plus a light ranking penalty.

The scSurvival reference mapping is documented in
`docs/fastde/scsurvival_reference_mapping.md`.
