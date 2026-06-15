# FastCopyKAT Project Plan

FastCopyKAT is a Python rewrite of the main CopyKAT execution contract. The first version keeps the public-facing artifacts compatible with CopyKAT while replacing several expensive R steps with vectorized NumPy/SciPy approximations.

## Compatibility Target

- Input: raw UMI counts in CopyKAT orientation, genes by cells. H5AD input is also accepted and converted to the same orientation.
- Required reference: gene coordinates with gene, chromosome, and position columns. CopyKAT-style `full.anno` works directly after export.
- Main outputs:
  - `<sample>_copykat_prediction.txt`
  - `<sample>_fastcopykat_cell_scores.tsv`
  - `<sample>_fastcopykat_manifest.json`

## Implemented Pipeline

1. Filter cells and genes using CopyKAT-style detection thresholds.
2. Transform counts with the same `log(sqrt(x) + sqrt(x + 1))` variance stabilizer.
3. Smooth per chromosome with a centered rolling filter.
4. Estimate a diploid baseline from known normal cells or from the flattest hierarchical cluster.
5. Segment smoothed relative expression with chromosome-local rolling differences.
6. Aggregate gene-level CNV profiles to 220 kb bins.
7. Predict diploid versus aneuploid cells from bin-level CNV amplitude.
8. Rescue local aneuploid calls with chromosome-level CNV burden.

## Simplifications

- The dynamic linear model smoother is approximated by a vectorized chromosome-local rolling mean.
- The MCMC segmentation is approximated by rolling-difference breakpoints and segment medians.
- Baseline selection keeps the CopyKAT cluster idea but uses robust CNV amplitude instead of mixture-model confidence tests.
- The default output is compact and skips the wide CNA TSV. Use `--output-mode copykat-tsv` only when the compatibility matrix is required.

These simplifications are intentional for the first benchmarkable version. They should be judged by agreement on test data and runtime, not by exact byte-level equality to the R internals.

## Benchmark Protocol

Use CopyKAT's own fixture for a first smoke comparison:

```bash
Rscript scripts/fastcopykat/export_copykat_fixture.R /tmp/copykat_src tmp/fastcopykat_fixture
PYTHONPATH=src .venv-cnmf-h2/bin/python -m fastcopykat run \
  --input tmp/fastcopykat_fixture/copykat_exp_rawdata.tsv \
  --gene-annotation tmp/fastcopykat_fixture/copykat_full_anno_hg20.tsv \
  --bins tmp/fastcopykat_fixture/copykat_DNA_hg20_bins.tsv \
  --output-dir tmp/fastcopykat_fixture/output \
  --sample-name test \
  --output-mode compact \
  --min-gene-per-cell 200
```

Compare `tmp/fastcopykat_fixture/output/test_copykat_prediction.txt` against `/tmp/copykat_src/test_output/test_copykat_prediction.txt`.

## Current Fixture Result

Dataset: CopyKAT package fixture `exp.rawdata`, 33,694 genes by 302 cells. Reference labels are from an actual local CopyKAT 1.2.3 run installed into `/tmp/Rlib-fastcopykat`.

| Engine | Wall time | Max RSS | Prediction counts |
| --- | ---: | ---: | --- |
| CopyKAT 1.2.3, `n.cores=1` | 106.72 s | 1.45 GB | 229 aneuploid, 73 diploid |
| FastCopyKAT Python, compact output | 6.22 s | 0.49 GB | 229 aneuploid, 73 diploid |

Observed speedup on this fixture: `17.2x` wall-clock. Peak memory decreased by about `66%`.

Prediction agreement against the local CopyKAT run:

| CopyKAT label | FastCopyKAT aneuploid | FastCopyKAT diploid |
| --- | ---: | ---: |
| aneuploid | 229 | 0 |
| diploid | 0 | 73 |

Overall agreement: `100.00%` across 302 overlapping cells. The compact output writes prediction, per-cell CNV scores, and manifest; the compatibility CNA TSV remains available with `--output-mode copykat-tsv`.

## Reproduction Notes

The R baseline was run after installing CopyKAT and its missing dependencies into a temporary library:

```bash
R CMD INSTALL -l /tmp/Rlib-fastcopykat /tmp/copykat_src
```

The Python run writes a manifest with per-stage timings. In the compact fixture run, internal compute time was 4.04 s; process wall time was 6.22 s including input parsing and compact output.

Use the comparison helper:

```bash
PYTHONPATH=src .venv-cnmf-h2/bin/python scripts/fastcopykat/compare_predictions.py \
  --reference tmp/fastcopykat_fixture/r_copykat/test_r_copykat_prediction.txt \
  --candidate tmp/fastcopykat_fixture/output_rescue_compact/test_copykat_prediction.txt \
  --output-json tmp/fastcopykat_fixture/fastcopykat_vs_copykat_prediction.json
```

## Efficiency Notes

The current acceleration comes from three places:

- Smoothing uses a vectorized chromosome-local rolling filter instead of per-cell dynamic linear model smoothing.
- Segmentation uses rolling-difference breakpoints and segment medians instead of MCMC sampling.
- Prediction avoids the transport/EMD-heavy branch and uses robust CNV amplitude plus chromosome-level burden rescue around the inferred diploid baseline.

The current bottleneck is no longer core compute on this fixture. It is mostly text I/O for the large CNA matrix. For larger data, the next targets are sparse/H5AD input streaming, Parquet/Zarr output, and optional compiled kernels for bin aggregation.
