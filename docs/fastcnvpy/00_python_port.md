# FastCNVpy Python Port

`fastcnvpy` is a Python port of the core fastCNV CNV calling path from
`must-bioinfo/fastCNV` 1.1.10. The implementation is intentionally conservative:
it mirrors the R matrix path first, then leaves performance work for later
optimization passes.

## Scope

Implemented:

- `CNVCalling`-compatible count normalization, gene selection, genomic windows,
  raw genomic scores, quantile trimming, and `cnv_fraction`.
- `CNVPerChromosomeArm`-style 46-arm summaries.
- `CNVclassification`-style arm calls with default peaks `-0.1, 0, 0.1`.
- Hierarchical CNV clustering and simple high-correlation cluster merging.
- H5AD, TSV, CSV, and Parquet inputs.
- Compact outputs for large runs, plus full TSV/Parquet matrix outputs for
  parity checks.

Not yet implemented:

- Seurat object mutation, Seurat plotting wrappers, ComplexHeatmap plots, tree
  plotting, Banksy spatial aggregation, Visium HD wrappers, and biomaRt online
  refresh logic.

The port uses the exported Ensembl v113 fastCNV `geneMetadata` table stored at
`tmp/fastcnvpy_reference/geneMetadata_ensembl113.tsv` for local tests.

## CLI

```bash
PYTHONPATH=src python -m fastcnvpy run \
  --input h5ad/canonical/quick_test/private_overall_sim_10samples_1000cells.h5ad \
  --gene-metadata tmp/fastcnvpy_reference/geneMetadata_ensembl113.tsv \
  --output-dir tmp/fastcnvpy_project_smoke \
  --sample-name private_quick \
  --output-mode compact \
  --no-clusters
```

Use `--output-mode tsv` or `--output-mode parquet` when exact matrix-level
inspection is needed. Compact mode writes only the manifest, per-cell metadata,
arm CNV matrix, and genomic window metadata.

For H5AD input, `--h5ad-mode dense` is the default and is the fastest mode.
`--h5ad-mode sparse` keeps the full expression matrix sparse until the selected
fastCNV genes are extracted; it reduces peak memory but is slower on the current
test H5AD.

For merged tumor cohorts, use the list-compatible pooled-reference interface:

```bash
PYTHONPATH=src python -m fastcnvpy run-pooled \
  --input merged_tumor.h5ad \
  --gene-metadata tmp/fastcnvpy_reference/geneMetadata_ensembl113.tsv \
  --output-dir fastcnvpy_pooled \
  --sample-key sample_id \
  --reference-var fastcnv_reference_pool \
  --reference-label normal_nonparenchymal,normal_parenchymal \
  --sample-name tumor_pooled \
  --h5ad-mode dense \
  --n-jobs 1 \
  --output-mode compact
```

This mirrors fastCNV's list input behavior: common genes and selected genes are
computed across samples, the scale factor and threshold quantiles are derived
from the pooled reference, and each sample is scored separately.

## Outputs

- `<sample>_fastcnvpy_manifest.json`
- `<sample>_fastcnvpy_cell_metadata.tsv`
- `<sample>_fastcnvpy_arm_cnv.tsv`
- `<sample>_fastcnvpy_genomic_windows.tsv`
- `<sample>_rawGenomicScores.tsv|parquet` when not using compact mode
- `<sample>_genomicScores.tsv|parquet` when not using compact mode

## R Parity Check

The standalone reference script
`scripts/fastcnvpy/run_r_cnvcalling_reference.R` reproduces the fastCNV
`CNVCalling` matrix path without requiring Seurat. On the smoke dataset:

- no-reference branch: raw and trimmed genomic scores matched the Python output
  with max absolute difference `8.88e-16`; `cnv_fraction` correlation was `1.0`.
- reference-label branch: raw and trimmed genomic scores matched with max
  absolute difference `8.88e-16`; `cnv_fraction` correlation was `1.0`.

One fastCNV/R-specific behavior is preserved deliberately: duplicated row names
in a matrix are resolved by character lookup to the first matching row. Pandas
normally returns all duplicated rows, so `fastcnvpy` uses first-match lookup when
calculating genomic window means.

## Performance

The first port prioritized exactness and was slower than the standalone R matrix
reference on medium inputs. The current implementation removes non-semantic
overhead while preserving identical Python outputs:

- no deep copy of the full counts matrix at pipeline start.
- numeric TSV inputs skip a second full-table `to_numeric` pass.
- average-expression and genomic-window score calculations use NumPy arrays
  instead of repeated DataFrame indexing.
- H5AD input uses an AnnData-specific path that avoids building a full
  gene-by-cell DataFrame. Dense mode computes expression means without copying a
  `cells x common_genes` submatrix.
- sparse H5AD mode is available for lower memory but is not the default speed
  path.
- pooled H5AD mode keeps the input merged for reference construction, then
  scores samples independently with optional threaded per-sample workers.
- duplicated gene names still follow R matrix first-match lookup.

Measured with clustering disabled:

| input | mode | Python before | Python after | R reference | peak RSS after |
| --- | --- | ---: | ---: | ---: | ---: |
| 30592 genes x 1000 cells TSV | full TSV matrices | 32.68s | 6.49s | 14.10s | 896 MB |
| 30592 genes x 10000 cells H5AD, dense mode | compact | 177.3s | 11.48s | not run | 4.13 GB |
| 30592 genes x 10000 cells H5AD, sparse mode | compact | not measured | 46.64s | not run | 2.43 GB |

For the 1000-cell TSV benchmark, optimized Python output was exactly identical
to the pre-optimization Python output and remained numerically equivalent to the
R reference (`raw` max absolute difference `9.99e-16`, `trimmed` max absolute
difference `9.44e-16`).

For the 10000-cell H5AD benchmark, dense AnnData output was exactly identical to
the previous dense DataFrame path for `cell_metadata` and `arm_cnv`. Sparse H5AD
mode was also exactly identical on those outputs, but was slower because sparse
column slicing dominated runtime on this dataset.

The remaining large memory cost in the fastest H5AD mode comes from densifying
the full expression matrix. On this benchmark, the low-memory sparse path cuts
peak RSS by about 41% but is about 4.1x slower than dense mode.
