# FastCopyKAT Scale Benchmark

This benchmark uses simulated scRNA count matrices built from CopyKAT hg20 gene symbols and genomic coordinates. The truth labels are known by construction: 30% diploid cells and 70% aneuploid cells with chromosome-scale gains and losses.

The goal is to test larger matrix behavior than the 302-cell CopyKAT package fixture while still measuring accuracy, runtime, and resource usage.

## Commands

FastCopyKAT 1k compact-output run:

```bash
/usr/bin/time -v env PYTHONPATH=src .venv-cnmf-h2/bin/python scripts/fastcopykat/benchmark_scale.py \
  --gene-annotation tmp/fastcopykat_fixture/copykat_full_anno_hg20.tsv \
  --bins tmp/fastcopykat_fixture/copykat_DNA_hg20_bins.tsv \
  --output-dir tmp/fastcopykat_scale_io \
  --sizes 1000 \
  --genes 12000 \
  --write-outputs \
  --output-mode compact \
  --save-inputs
```

CopyKAT 1k reference run:

```bash
/usr/bin/time -v Rscript scripts/fastcopykat/run_copykat_reference_synthetic.R \
  tmp/fastcopykat_scale_io/synthetic_1000/synthetic_1000_counts.tsv \
  tmp/fastcopykat_scale_io/r_copykat_synthetic_1000 \
  synthetic_1000_r \
  1
```

FastCopyKAT 3k/5k compact-output run:

```bash
/usr/bin/time -v env PYTHONPATH=src .venv-cnmf-h2/bin/python scripts/fastcopykat/benchmark_scale.py \
  --gene-annotation tmp/fastcopykat_fixture/copykat_full_anno_hg20.tsv \
  --bins tmp/fastcopykat_fixture/copykat_DNA_hg20_bins.tsv \
  --output-dir tmp/fastcopykat_rescue_scale_compact \
  --sizes 3000 5000 \
  --genes 12000 \
  --write-outputs \
  --output-mode compact
```

## Results

| Dataset | Engine | Output mode | Accuracy | Balanced accuracy | Compute | Write | Wall | Max RSS | Output size |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1k cells x 12k genes | FastCopyKAT rescue | compact | 100.00% | 100.00% | 11.31 s | 0.03 s | 14.54 s process / 11.34 s pipeline | 0.94 GB | 128 KB |
| 1k cells x 12k genes | CopyKAT 1.2.3 | raw + final CNA TSV + heatmap | 100.00% | 100.00% | n/a | n/a | 408.08 s | 3.78 GB | 665 MB |
| 3k cells x 12k genes | FastCopyKAT rescue | compact | 100.00% | 100.00% | 63.23 s | 0.04 s | n/a / 63.27 s pipeline | 2.60 GB | 380 KB |
| 5k cells x 12k genes | FastCopyKAT rescue | compact | 100.00% | 100.00% | 89.05 s | 0.06 s | n/a / 89.12 s pipeline | 4.24 GB | 628 KB |

For the 1k compact-output comparison, FastCopyKAT was about `28.1x` faster end-to-end than CopyKAT and used about `75%` less peak memory. With chromosome-level rescue, FastCopyKAT recovered the planted labels on this synthetic signal.

## 1k Confusion Against Truth

FastCopyKAT:

| Truth | Predicted aneuploid | Predicted diploid |
| --- | ---: | ---: |
| aneuploid | 700 | 0 |
| diploid | 0 | 300 |

CopyKAT:

| Truth | Predicted aneuploid | Predicted diploid |
| --- | ---: | ---: |
| aneuploid | 700 | 0 |
| diploid | 0 | 300 |

Engine agreement on the 1k benchmark was `100.00%`.

## Bottlenecks

At larger scale the dominant costs are now mostly compute-side when compact output is used.

- Compact 5k: baseline selection dominates. The 5k run spent 74.05 s of 89.04 s internal time in hierarchical baseline estimation.
- Compact output removes the wide TSV bottleneck. The 3k run spent 0.04 s writing compact outputs, versus 166.22 s for the previous 727 MB final CNA TSV.
- Memory still scales with dense gene-by-cell and bin-by-cell matrices. The 5k compact run peaked at 4.24 GB.

## Interpretation

FastCopyKAT is clearly faster and lighter than CopyKAT on these benchmarks. The chromosome-level rescue rule fixed the previous aneuploid sensitivity loss on the planted synthetic data and also matched the 302-cell CopyKAT fixture labels exactly. This does not prove biological equivalence on all tumors; it shows the earlier loss was caused by over-dilution in whole-genome median CNV scoring.

For speed, the next high-impact changes are:

- Replace hierarchical baseline selection with MiniBatchKMeans or approximate nearest-neighbor clustering for large `n_cells`.
- Keep compact output as the default; store final CNA as Parquet/Zarr/H5AD only when users need the matrix.
- Stream bin aggregation and output in chunks to reduce peak memory.
