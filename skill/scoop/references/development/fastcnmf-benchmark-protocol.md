# FastCNMF Benchmark Protocol

## Baselines

FastCNMF benchmarks compare three lanes:

1. `cnmf-reference`: upstream `cnmf==1.7.1` plus Harmony-compatible
   preprocessing.
2. `fastcnmf-exact-compatible`: same math contract, improved orchestration and
   caching.
3. `fastcnmf-accelerated`: GPU, float32, dynamic scheduling, or approximation
   enabled as configured.

## Dataset Tiers

| Tier | Dataset | Purpose |
| --- | --- | --- |
| S0 | 3 selected GBM Visium samples | smoke test and CI-scale runtime |
| S1 | all available GBM Visium samples | spatial medium benchmark |
| S2 | 72k public + 72k internal scRNA fixtures | medium single-cell benchmark |
| S3 | 200k-500k sampled cells/spots | large benchmark |
| S4 | 1M cells/spots sketch/refit | stress test for approximate mode |

## Fixed Parameters

Unless a benchmark explicitly studies parameters:

- Harmony batch key: `sample_id`
- HVG count: `3000`
- k values: `6, 8, 10, 12`
- NMF replicates per k: `20` for medium/large, `8` for smoke
- cNMF/FastCNMF max NMF iterations: `50` by default for smoke and production-scale tests
- worker counts: `1, 4, 8`, plus GPU-specific device counts where available

## Runtime and Resource Metrics

Collect for each stage:

- wall-clock seconds
- user/system CPU seconds
- CPU utilization
- max RSS
- file-system input/output counts
- output directory size
- GPU peak allocated/reserved VRAM
- GPU utilization and OOM retries, if applicable

## Quality Metrics

Compare FastCNMF outputs to cNMF baseline after matching programs by Hungarian
assignment on spectra cosine similarity.

- spectra cosine similarity
- usage Pearson correlation
- usage Spearman correlation
- top gene Jaccard index
- reconstruction error and relative delta
- k-selection stability ranking
- sample/batch leakage in usage, using `sample_id` predictability
- optional spatial smoothness metrics for spatial inputs

## GPU OOM Test

For GPU backends, run a matrix-size sweep:

| Sweep Axis | Values |
| --- | --- |
| observations | `10k, 50k, 100k, 250k, 500k` |
| genes | `2k, 3k, 5k` |
| k | `8, 16, 32` |
| concurrent replicates | `1, 2, 4, 8` |
| dtype | `float32`, optional `float16/bfloat16` |

An OOM test passes only if FastCNMF:

- catches OOM at task level
- records the failed estimate and actual peak memory
- retries with a smaller batch/chunk or CPU fallback
- does not corrupt completed replicate outputs

## Success Criteria

### Exact-Compatible

- factorize speedup `>= 1.5x`
- end-to-end speedup `>= 1.2x`
- spectra cosine `>= 0.995`
- usage correlation `>= 0.99`

### GPU Fast-Compatible

- factorize speedup `>= 5x`
- end-to-end speedup `>= 2x`
- spectra cosine `>= 0.98`
- usage correlation `>= 0.95`
- reconstruction error delta `<= 2%`

### Large Approximate

- completes datasets where cNMF baseline is infeasible under resource budget
- spectra cosine on validation subset `>= 0.95`
- reconstruction error delta on validation subset `<= 5%`

## Required Report Artifacts

Each benchmark run writes:

- `benchmark_summary.json`
- `benchmark_report.md`
- `resource_timeseries.tsv`
- `environment_freeze.txt`
- `task_manifest.json`
- `quality_metrics.json`
- cNMF-compatible exported usages/spectra files
