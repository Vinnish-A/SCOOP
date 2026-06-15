# FastCopyKAT Project Sample Benchmark

This benchmark splits the project quick-test H5AD files by `sample_id`, samples 300 cells from each selected sample, and runs each sample independently through FastCopyKAT compact output and R CopyKAT 1.2.3.

Datasets:

- `public`: `h5ad/canonical/quick_test/public_O_GSE154795_10samples_1000cells.h5ad`
- `private`: `h5ad/canonical/quick_test/private_overall_sim_10samples_1000cells.h5ad`

Command:

```bash
PYTHONPATH=src .venv-cnmf-h2/bin/python scripts/fastcopykat/benchmark_project_samples.py \
  --samples-per-dataset 2 \
  --cells-per-sample 300 \
  --output-dir tmp/fastcopykat_project_sample_benchmark_300
```

The full generated report is at `tmp/fastcopykat_project_sample_benchmark_300/summary.md`.

## Results

| Dataset | Sample | Cells | Fast wall | CopyKAT wall | Speedup | Fast RSS | CopyKAT RSS | Defined agreement | Fast output | CopyKAT output |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| private | A_4843092 | 300 | 7.20 s | 521.13 s | 72.38x | 0.42 GB | 1.22 GB | 75.93% | 0.05 MB | 175.70 MB |
| private | A_4888223 | 300 | 7.84 s | 166.07 s | 21.18x | 0.41 GB | 1.22 GB | 77.66% | 0.04 MB | 162.75 MB |
| public | ndGBM-10 | 300 | 4.17 s | 88.00 s | 21.10x | 0.39 GB | 0.87 GB | 83.33% | 0.04 MB | 74.04 MB |
| public | rGBM-01-B | 300 | 10.09 s | 167.63 s | 16.61x | 0.39 GB | 1.22 GB | 77.89% | 0.06 MB | 156.17 MB |

Mean:

- FastCopyKAT wall time: `7.33 s`
- CopyKAT wall time: `235.71 s`
- Mean speedup: `32.82x`
- FastCopyKAT peak RSS: `0.40 GB`
- CopyKAT peak RSS: `1.13 GB`
- Defined normalized agreement: `78.70%`
- CopyKAT output was about `2,980x` larger than FastCopyKAT compact output on average.

## Agreement Definition

The quick-test samples do not include biological CNV truth labels, so this benchmark measures agreement with R CopyKAT, not biological accuracy.

R CopyKAT sometimes returns labels such as `c1:diploid:low.conf`, `c2:aneuploid:low.conf`, and `not.defined`. For the main agreement number:

- labels containing `diploid` are mapped to `diploid`
- labels containing `aneuploid` are mapped to `aneuploid`
- `not.defined` is excluded

Cells can also be filtered by either method, so overlap can be smaller than 300.

## Per-Sample Confusion

Rows are normalized CopyKAT labels and columns are FastCopyKAT labels.

`private/A_4843092`, 295 defined overlapping cells:

| CopyKAT | Fast aneuploid | Fast diploid |
| --- | ---: | ---: |
| aneuploid | 218 | 35 |
| diploid | 36 | 6 |

`private/A_4888223`, 273 defined overlapping cells:

| CopyKAT | Fast aneuploid | Fast diploid |
| --- | ---: | ---: |
| aneuploid | 170 | 56 |
| diploid | 5 | 42 |

`public/ndGBM-10`, 138 defined overlapping cells:

| CopyKAT | Fast aneuploid | Fast diploid |
| --- | ---: | ---: |
| aneuploid | 31 | 4 |
| diploid | 19 | 84 |

`public/rGBM-01-B`, 285 defined overlapping cells:

| CopyKAT | Fast aneuploid | Fast diploid |
| --- | ---: | ---: |
| aneuploid | 187 | 62 |
| diploid | 1 | 35 |

## Interpretation

FastCopyKAT is consistently much faster and lighter on per-sample project test subsets, mainly because it avoids CopyKAT's per-cell MCMC-heavy segmentation and wide text outputs. However, project-sample agreement is lower than on the CopyKAT package fixture and planted synthetic tests.

The disagreement pattern is sample-dependent:

- Some samples show FastCopyKAT calling fewer aneuploid cells than CopyKAT.
- One private sample shows FastCopyKAT calling more aneuploid cells among CopyKAT diploid cells.
- CopyKAT emits many `not.defined` cells for `public/ndGBM-10`, so direct string agreement is not meaningful for that sample.

The next accuracy improvement should use real project-sample calibration, not only planted synthetic data. A practical next step is to tune `chromosome_rescue_mad_multiplier` and add a low-confidence output state instead of forcing every retained cell into diploid/aneuploid.

