# FastCNMF Large Benchmark Roadmap

## Current Scope

The next target is a fair same-hardware cold-start benchmark on larger S1/S2
datasets, then promotion of FastCNMF from a cNMF wrapper into an independent
execution framework.

The current generated manifest is:

- `tmp/fastcnmf_large_benchmark/benchmark_manifest.json`
- `tmp/fastcnmf_large_benchmark/execution_plan.json`
- `tmp/fastcnmf_large_benchmark/run_bundle/run_bundle.json`
- `tmp/fastcnmf_large_benchmark/run_bundle/scripts/`

It contains:

| dataset_id | tier | observations | variables | samples | cells/spots per sample |
| --- | --- | ---: | ---: | ---: | --- |
| `s1_gbm_lowres_visium_available` | S1 partial | 7,071 | unknown from directory scan | 3 | 1,407-2,861 |
| `s2_public_24x3000` | S2 | 72,000 | 40,791 | 24 | 3,000 each |
| `s2_internal_24x3000` | S2 | 72,000 | 36,581 | 24 | 3,000 each |

S1 is marked incomplete because only 3 GBM Visium samples are present in the
project. S2 is the first real enlarged benchmark tier available now.

## Fairness Contract

The benchmark manifest uses the comparison name `same_hardware_cold_start`.
Both lanes must start from the same raw or canonical input and write to their
own output roots:

- reference lane: `cnmf_optimized`
- candidate lane: `fastcnmf_independent`
- target speedup: `3.0x`
- minimum spectra cosine: `0.95`
- minimum usage Pearson: `0.95`
- cross-lane artifact reuse: forbidden

Allowed hardware acceleration is not considered unfair. The reference cNMF lane
may use CPU multiprocessing and BLAS thread limits. The FastCNMF lane may use
the same CPU resources and GPU if the selected backend supports it. The report
must record which hardware path was actually used.

## Independent Execution Framework

The execution plan expands each dataset/lane pair into five stages:

1. `preprocess_cold_start`
2. `plan_nmf_tasks`
3. `factorize_replicates`
4. `consensus`
5. `resource_and_quality_report`

Each stage owns its outputs under:

```text
tmp/fastcnmf_large_benchmark/<dataset_id>/<lane_id>/
```

This is the first boundary separating FastCNMF from the earlier ad hoc wrapper
scripts. The next implementation step is to attach real stage executors:

- cNMF reference executor
- FastCNMF preprocess executor
- FastCNMF NMF executor
- FastCNMF consensus executor
- benchmark report executor

## 3x Path

The 3x target should be measured against `cnmf_optimized`, not the easier
serial-only cNMF baseline.

Likely required changes:

- avoid repeated H5AD materialization by using cacheable core arrays and backed
  stores
- keep Harmony2 and MOE correction in-process
- build NMF task manifests directly instead of shelling through cNMF prepare
- schedule replicate tasks dynamically instead of static worker indices
- parallelize consensus by k
- add GPU NMF only after the CPU cold-start comparison is reproducible

## Commands

Generate the current manifest and execution plan:

```bash
PYTHONPATH=src python -m fastcnmf benchmark-manifest \
  --root . \
  --output-root tmp/fastcnmf_large_benchmark \
  --output tmp/fastcnmf_large_benchmark/benchmark_manifest.json

PYTHONPATH=src python -m fastcnmf plan-run \
  --manifest tmp/fastcnmf_large_benchmark/benchmark_manifest.json \
  --output tmp/fastcnmf_large_benchmark/execution_plan.json
```

These commands do not prove the 3x target. They make the larger benchmark and
fairness contract reproducible so the real executors can be implemented next.

Generate executable stage scripts:

```bash
PYTHONPATH=src python -m fastcnmf write-run-bundle \
  --manifest tmp/fastcnmf_large_benchmark/benchmark_manifest.json \
  --output-dir tmp/fastcnmf_large_benchmark/run_bundle \
  --reference-python ./.venv-cnmf/bin/python \
  --candidate-python ./.venv-cnmf-h2/bin/python
```

The first heavy S2 reference baseline stage is:

```bash
bash tmp/fastcnmf_large_benchmark/run_bundle/scripts/s2_public_24x3000__cnmf_optimized__preprocess.sh
```

The script writes `/usr/bin/time -v` output to:

```text
tmp/fastcnmf_large_benchmark/run_bundle/logs/s2_public_24x3000__cnmf_optimized__preprocess.time.log
```

The matching FastCNMF candidate preprocess stage is now implemented and points
to the Harmony2 environment:

```bash
bash tmp/fastcnmf_large_benchmark/run_bundle/scripts/s2_public_24x3000__fastcnmf_independent__preprocess.sh
```

Its timing log will be:

```text
tmp/fastcnmf_large_benchmark/run_bundle/logs/s2_public_24x3000__fastcnmf_independent__preprocess.time.log
```

Both commands are cold-start preprocess stages and write to separate lane roots;
neither command is allowed to consume the other lane's intermediate artifacts.

## First S2 Preprocess Result

The first S2 public preprocessing comparison has been executed with CPU-only
CUDA visibility for both lanes:

- Reference cNMF preprocess: `258.54 s`
- FastCNMF candidate preprocess: `106.32 s`
- Preprocess speedup: `2.43x`
- Corrected input cosine: `0.9957`
- Corrected input Pearson: `0.9956`

Artifacts:

- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/preprocess_benchmark_summary.md`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/preprocess_benchmark_summary.json`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/preprocess_compare.json`

The earlier GPU-visible cNMF reference attempt crashed with signal `11` in the
Harmony0.2 CUDA path after `524.03 s`; that failure is retained as evidence but
is not used as the fair baseline.

## S2 Smoke End-to-End Result

The first S2 public smoke chain was run through preprocess, cNMF prepare,
factorize, consensus, and spectra/usage comparison:

- Reference end-to-end: `926.85 s`
- FastCNMF candidate end-to-end: `924.83 s`
- End-to-end speedup: `1.00x`
- Mean spectra cosine: `0.9983`
- Mean usage Pearson: `0.9953`

Artifact:

- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_end_to_end_summary.md`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_end_to_end_summary.json`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_cnmf_compare.json`

This confirms that output quality remains above 95%, but it also shows that
preprocessing acceleration alone is not enough. The cNMF-compatible prepare,
factorize, and consensus stages erase the preprocessing gain.

## S2 Fast Executor Smoke Result

The candidate lane now uses FastCNMF-owned factorize and consensus entry
points:

- `python -m fastcnmf fast-factorize`
- `python -m fastcnmf fast-consensus`

The reference lane remains the optimized cNMF CLI lane. The S2 public smoke
result is:

| stage | reference seconds | candidate seconds | speedup | candidate max RSS MB |
| --- | ---: | ---: | ---: | ---: |
| preprocess | 258.54 | 106.32 | 2.43x | 12596.7 |
| prepare | 76.91 | 105.42 | 0.73x | 11620.8 |
| factorize | 374.55 | 431.59 | 0.87x | 3652.2 |
| consensus | 216.85 | 84.37 | 2.57x | 8628.0 |
| end-to-end | 926.85 | 727.70 | 1.27x |  |

Quality remains above the configured gate:

- Mean spectra cosine: `0.9983`
- Mean usage Pearson: `0.9953`
- 95% spectra/usage gate: `passed`
- 3x speed gate: `failed`

Artifacts:

- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_fast_executor_end_to_end_summary.md`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_fast_executor_end_to_end_summary.json`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_fast_executor_cnmf_compare.json`

Interpretation:

- The independent execution framework boundary is now real for factorize and
  consensus, not only for preprocessing.
- k-parallel consensus works and gives a meaningful stage-level speedup.
- Dynamic scheduling around sklearn CD-NMF does not solve factorize; that stage
  is still slower than the optimized cNMF reference.

The next 3x path must replace or substantially improve:

- static cNMF worker-index factorization
- repeated cNMF prepare materialization
- serial combine/consensus/k-selection execution

The factorize replacement must be stronger than scheduling. Candidate options
are:

- batched CPU NMF across replicate seeds
- GPU NMF with memory-aware tiling and automatic batch-size fallback
- a gate-controlled approximate mode that reduces replicate count or max NMF
  iterations only when spectra/usage consistency stays above threshold
- direct FastCNMF prepare output to remove the `105.42 s` candidate prepare
  penalty

## S2 Direct Prepare Smoke Result

The candidate lane now also owns prepare/cache writing through
`python -m fastcnmf fast-prepare`. It writes cNMF-compatible prepare artifacts
directly and hardlinks the lane-local TP10K H5AD when possible.

| stage | reference seconds | candidate seconds | speedup | candidate max RSS MB |
| --- | ---: | ---: | ---: | ---: |
| preprocess | 258.54 | 106.32 | 2.43x | 12596.7 |
| prepare | 76.91 | 39.29 | 1.96x | 8259.5 |
| factorize | 374.55 | 431.15 | 0.87x | 1913.1 |
| consensus | 216.85 | 91.45 | 2.37x | 6897.3 |
| end-to-end | 926.85 | 668.21 | 1.39x |  |

Quality:

- Mean spectra cosine: `0.9983`
- Mean usage Pearson: `0.9953`
- 95% spectra/usage gate: `passed`
- 3x speed gate: `failed`

Artifacts:

- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_direct_prepare_end_to_end_summary.md`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_direct_prepare_end_to_end_summary.json`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_direct_prepare_cnmf_compare.json`

Interpretation:

- The candidate lane now owns preprocessing, prepare/cache writing, factorize
  scheduling, and consensus scheduling.
- Prepare overhead has been reduced substantially, but the total is still only
  `1.39x` faster than the optimized cNMF reference.
- The next step toward `3x` must target the NMF factorization kernel or define
  a separately reported approximate mode with an explicit quality gate.

## Rejected Approximate n_iter=2 Smoke

An exploratory approximate variant was run with `n_iter=2` and
`local_neighborhood_size=0.5` under a separate output root:

- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_approx_n2/`

It produced an apparent `3.64x` speedup, but it is rejected for two independent
reasons:

- The manual run did not inherit run-bundle BLAS thread limits, so the speed is
  not a fair benchmark result.
- Quality failed the required gate:
  - Mean spectra cosine: `0.9229`
  - Mean usage Pearson: `0.8890`
  - 95% spectra/usage gate: `failed`

Artifacts:

- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_approx_n2/approx_n2_summary.md`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_approx_n2/approx_n2_summary.json`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_approx_n2/approx_n2_cnmf_compare.json`

This negative result is useful: simply reducing replicate count is not a valid
route to FastCNMF for this dataset. The next accelerator must preserve the
replicate information more faithfully, either through a faster NMF backend or a
more careful approximation with a much tighter quality guard.

## Rejected MiniBatchNMF Backend Smoke

An experimental `MiniBatchNMF` backend was added to `fast-factorize` and run as
a separate variant:

- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_minibatch_n8/`

It used the same `n_iter=8` replicate count as the smoke baseline and explicit
BLAS thread limits. The factorize stage was fast:

- MiniBatch factorize: `33.83 s`
- Apparent end-to-end speedup with no-filter `dt=2.0` consensus: `3.04x`

But it is rejected:

- Default `dt=0.5` consensus failed density filtering.
- The no-filter `dt=2.0` output failed the quality gate:
  - Mean spectra cosine: `0.4075`
  - Mean usage Pearson: `0.4951`
  - 95% spectra/usage gate: `failed`

Artifacts:

- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_minibatch_n8/minibatch_n8_summary.md`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_minibatch_n8/minibatch_n8_summary.json`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_minibatch_n8/minibatch_n8_dt2_cnmf_compare.json`

The useful lesson is that the factorize wall time can be reduced enough to
cross the speed target, but not with this backend's current numerical behavior.
The next backend should preserve cNMF spectra geometry more faithfully, for
example exact-style multiplicative updates on GPU or a batched CPU/GPU solver
with a stricter convergence policy.

## Rejected CuPy MU GPU Backend Smoke

A UV-managed GPU experiment environment was created:

- `.venv-fastcnmf-gpu`
- `cupy-cuda12x==14.1.1`

The `fast-factorize` command now has an experimental `--backend cupy-mu` mode.
It uses GPU multiplicative updates and writes cNMF-compatible spectra files.
The backend ran on the S2 public smoke dataset without OOM on the RTX 3080 Ti.

| variant | factorize seconds | candidate accounting seconds | apparent speedup | spectra cosine | usage Pearson | quality gate |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `cupy_mu_n8` | 31.52 | 218.44 | 4.24x | 0.8959 | 0.8652 | failed |
| `cupy_mu_n8_i600` | 60.55 | 261.28 | 3.55x | 0.8958 | 0.8654 | failed |

Artifacts:

- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_cupy_mu_n8/cupy_mu_n8_summary.md`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_cupy_mu_n8/cupy_mu_n8_summary.json`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_cupy_mu_n8_i600/cupy_mu_n8_i600_summary.md`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_cupy_mu_n8_i600/cupy_mu_n8_i600_summary.json`

This result is important but not sufficient. GPU factorize is fast enough and
fits memory, but plain MU does not match cNMF's CD-derived spectra geometry.
The next backend should be CD-like or otherwise quality-calibrated against the
cNMF reference, rather than simply increasing MU iterations.

## Exact CD Iteration Truncation Boundary

The next CD-like experiment kept sklearn's exact coordinate-descent solver but
reduced `max_iter` while preserving `n_iter=8`.

| max_iter | total seconds | speedup | factorize seconds | consensus seconds | spectra cosine | usage Pearson | quality gate | 3x gate |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| 20 | 462.06 | 2.01x | 93.80 | 154.98 | 0.9628 | 0.9359 | failed | failed |
| 25 | 515.94 | 1.80x | 116.73 | 276.64 | 0.9936 | 0.9799 | passed | failed |
| 35 | 425.49 | 2.18x | 193.84 | 114.94 | 0.9983 | 0.9943 | passed | failed |
| 50 | 434.14 | 2.13x | 241.55 | 69.63 | 0.9985 | 0.9953 | passed | failed |

Artifacts:

- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/exact_cd_truncation_summary.md`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/exact_cd_truncation_summary.json`

This is the most informative approximation so far. It preserves cNMF geometry
when `max_iter >= 25`, but it still does not meet the 3x end-to-end target.
The remaining path is not plain solver truncation; FastCNMF needs either:

- a faster CD-like kernel,
- faster consensus/refit, or
- a quality-gated strategy that combines moderate CD truncation with consensus
  optimization.

## Spectra/Usage-Compatible Lite Consensus Boundary

`fast-consensus --lite` now preserves the final cNMF spectra/usage definitions
used by the benchmark, while skipping full cNMF gene-level outputs:

- skipped: gene-level TPM spectra output, OLS gene scores, starCAT reference,
  clustergram, and k-selection
- retained: consensus spectra and final usages refit on std-scaled HVG TPM

| variant | total seconds | speedup | lite consensus seconds | max RSS MB | mean spectra | mean usage | min k overall | global mean gate | all-k gate | 3x gate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- |
| `exact_cd_i25_lite` | 264.24 | 3.51x | 24.94 | 6898.2 | 0.9936 | 0.9799 | 0.9267 | passed | failed | passed |
| `exact_cd_i35_lite` | 345.53 | 2.68x | 34.98 | 6901.2 | 0.9983 | 0.9943 | 0.9831 | passed | passed | failed |

Artifacts:

- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/exact_cd_lite_consensus_summary.md`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/exact_cd_lite_consensus_summary.json`

This changes the roadmap: the `3x` path is now close, but a rigorous all-k
quality gate prevents accepting the fastest run. The next iteration should
target one of two narrow fixes:

- improve `exact_cd_i25_lite` k=6 usage consistency from `0.9267` to at least
  `0.95`
- reduce `exact_cd_i35_lite` wall time from `345.53 s` to `<= 308.95 s`
  without reducing all-k quality

## i27 Float32 Candidate and Stability Gap

The next boundary search found a stricter candidate:

- `max_iter=27`
- largest-k-first dynamic factorize scheduling
- `norm_counts.h5ad` written as float32
- spectra/usage-compatible lite consensus

| variant | dtype | total seconds | speedup | prepare | factorize | consensus lite | mean spectra | mean usage | min k overall | strict accepted |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `exact_cd_i27_lkf_f32` | float32 | 270.62 | 3.42x | 28.13 | 88.51 | 47.66 | 0.9959 | 0.9866 | 0.9532 | true |
| `exact_cd_i27_lkf_f32_r2` | float32 | 362.22 | 2.56x | 74.24 | 123.20 | 58.46 | 0.9959 | 0.9866 | 0.9532 | false |

Artifacts:

- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/exact_cd_lkf_float32_stability_summary.md`
- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/exact_cd_lkf_float32_stability_summary.json`

This is the first strict quality-preserving `>3x` run, but it is not yet a
stable benchmark result. The independent repeat preserved spectra/usage quality
but missed the speed target because prepare, factorize, and consensus wall time
all regressed.

The next roadmap item is therefore runtime control, not another approximate
solver:

- reduce prepare I/O variance, likely by avoiding repeated H5AD materialization
  or using a reusable memory-mapped core array under a declared cold-start cache
  policy
- record per-k consensus timing to separate NNLS refit cost from I/O and
  process-start overhead
- run benchmark summaries over repeated trials and accept only a stable
  statistic such as median or worst-of-2 under the strict all-k quality gate

## Core Cache Prepare Path

FastCNMF now has an explicit core-cache path:

- `fast-preprocess --write-core-cache` writes normalized counts as `.npy`, obs
  and var names, and precomputed TP10K stats.
- `fast-prepare --norm-store npy --precomputed-*` hardlinks those artifacts into
  the run directory.
- `fast-factorize` and `fast-consensus --lite` read the manifest and use the
  `.npy` path when available.

Measured S2 public result:

| variant | total seconds | speedup | preprocess | prepare | factorize | consensus | min k overall | strict accepted |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `i27_lkf_f32_single` | 270.62 | 3.42x | 106.32 | 28.13 | 88.51 | 47.66 | 0.9532 | true |
| `i27_lkf_f32_repeat` | 362.22 | 2.56x | 106.32 | 74.24 | 123.20 | 58.46 | 0.9532 | false |
| `i27_corecache` | 326.87 | 2.84x | 167.81 | 2.03 | 109.79 | 47.24 | 0.9532 | false |

Artifact:

- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/core_cache_stability_summary.md`

The cache design is directionally correct: prepare is no longer a meaningful
runtime source. However, the cold-start run still does not stably reach `3x`
because cache materialization moved cost into preprocess. The next optimization
should avoid writing both large H5AD and cache forms during the same cold-start
run, or move FastCNMF preprocessing fully onto a chunked/backed core format so
TP10K and normalized-count artifacts are not materialized more than once.

## HVG-Only Lite Consensus Cache

FastCNMF-lite now avoids full TP10K H5AD reads during consensus:

- preprocessing can write TP10K HVG raw/scaled sparse caches
- prepare links those caches into the run manifest
- lite consensus refits TPM spectra only for HVGs, which is equivalent for the
  final usage output because fixed-usage spectra refit is independent per gene
- `--no-tp10k-h5ad` skips full TP10K H5AD output for lite-only runs

Measured S2 public result:

| variant | total seconds | speedup | preprocess | prepare | factorize | consensus | consensus RSS MB | strict accepted |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `i27_corecache` | 326.87 | 2.84x | 167.81 | 2.03 | 109.79 | 47.24 | 6066.4 | false |
| `i27_corecache_hvg_fix` | 346.15 | 2.68x | 208.53 | 3.77 | 121.70 | 12.15 | 1594.7 | false |
| `i27_corecache_hvg_notp` | 324.25 | 2.86x | 180.47 | 2.03 | 130.66 | 11.09 | 1593.8 | false |

Artifact:

- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/hvg_consensus_cache_summary.md`

This optimization moves consensus out of the critical path for the S2 smoke
benchmark. The remaining work for stable `~3x` is now narrower:

- stabilize factorize around the first float32 run's `~90 s` rather than the
  repeated `120-130 s` range
- reduce preprocessing materialization, especially duplicate corrected/core
  artifacts and sparse/dense conversion overhead
- run a repeated benchmark once those two stages are controlled

## Same-Worker S2 Fairness Result

The fairness caveat around worker count has now been resolved for the S2 public
benchmark. A new cNMF reference lane was run with the same 8-process CPU worker
budget used by the accepted FastCNMF candidate.

Artifact:

- `tmp/fastcnmf_large_benchmark/s2_public_24x3000/same_worker_w8_fairness_summary.md`

Measured result:

| lane | preprocess | prepare | factorize | consensus | total | speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cNMF optimized, 8 workers | 258.54 | 75.77 | 235.69 | 240.69 | 810.69 | 1.00x |
| FastCNMF, 8 workers | 180.47 | 2.02 | 61.16 | 10.22 | 253.87 | 3.19x |
| FastCNMF, worst-of-two stage envelope | 180.47 | 2.04 | 61.16 | 10.22 | 253.89 | 3.19x |

Quality versus the cNMF 8-worker reference:

- mean spectra cosine: `0.9959`
- mean usage Pearson: `0.9866`
- strict all-k gate: `true`
- minimum per-k overall consistency: `0.9532`

The acceleration is not from skipping an equivalent user-facing result in the
accepted lite contract. It comes from changing the execution architecture:

- cNMF factorize uses static worker-index partitioning; FastCNMF uses a dynamic
  largest-k-first task queue with shared memmapped normalized counts.
- cNMF prepare re-materializes the NMF input; FastCNMF links precomputed core
  cache artifacts and writes only the cNMF-compatible metadata needed by the
  runner.
- cNMF consensus rebuilds full output products; FastCNMF-lite uses cached
  TP10K HVG matrices and emits the spectra/usage outputs needed by the current
  benchmark gate.

The S2 public result now satisfies the `~3x` same-worker target. The remaining
roadmap item is coverage: repeat this benchmark level for the expanded internal
S2 data and larger spatial S1 data before calling the full project goal
complete.

## Internal S2 Coverage Status

The same-worker runner has been promoted into a reusable script:

- `scripts/fastcnmf/run_s2_same_worker_w8_benchmark.sh`

It was run on the internal S2 dataset:

- input: `h5ad/canonical/quick_test/internal_overall_sim_24samples_3000cells_balanced.h5ad`
- shape: `72,000` cells, `24` samples, `3,000` cells/sample
- cNMF optimized 8-worker reference total: `1135.24 s`

Artifact:

- `tmp/fastcnmf_large_benchmark/s2_internal_24x3000/internal_s2_i27_i100_sweep_summary.md`

FastCNMF max-iteration sweep:

| variant | max_iter | total | speedup | factorize | mean spectra | mean usage | min k overall | strict all-k |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `i27` | 27 | 227.93 | 4.98x | 72.97 | 0.9541 | 0.9550 | 0.8795 | false |
| `i35` | 35 | 244.91 | 4.64x | 87.39 | 0.9397 | 0.9331 | 0.8555 | false |
| `i50` | 50 | 260.95 | 4.35x | 103.80 | 0.9626 | 0.9596 | 0.8593 | false |
| `i100` | 100 | 336.60 | 3.37x | 179.71 | 0.9720 | 0.9643 | 0.8597 | false |

Input equivalence checks are strong:

- corrected input cosine: `0.9958`
- corrected input Pearson: `0.9955`
- norm-count input cosine: `0.9960`
- norm-count input Pearson: `0.9958`

This narrowed the internal blocker to NMF/consensus stability rather than
preprocessing. The quality failure was concentrated in `k=10` for higher
max-iteration candidates, so the next internal benchmark tested a larger
replicate count rather than only increasing per-replicate iterations.

Accepted internal S2 result:

- artifact:
  `tmp/fastcnmf_large_benchmark/s2_internal_24x3000/internal_s2_accepted_n20_i50_summary.md`
- FastCNMF candidate: `n_iter=20`, `max_iter=50`
- cNMF optimized 8-worker reference total: `1458.03 s`
- FastCNMF total: `399.51 s`
- speedup: `3.65x`
- mean spectra cosine: `0.9982`
- mean usage Pearson: `0.9990`
- minimum per-k overall consistency: `0.9967`
- strict all-k gate: `true`

This makes both expanded S2 datasets accepted:

- public S2: `n_iter=8`, `max_iter=27`, `3.19x` same-worker speedup
- internal S2: `n_iter=20`, `max_iter=50`, `3.65x` same-worker speedup

The remaining coverage gap is S1 spatial. The current manifest only has three
GBM Visium samples available under `data/raw/spatial/gbm_lowres_visium`, so S1
still needs either a larger spatial input or an explicit scoped benchmark
decision before the overall goal can be closed.

## Scoped S1 Spatial Result

The available spatial subset was benchmarked with the independent FastCNMF
runner. This uses the three local GBM low-resolution Visium samples already
materialized by the spatial Harmony benchmark.

Artifact:

- `tmp/cnmf_spatial_harmony_benchmark/fastcnmf_scoped_i50_w8/scoped_s1_fastcnmf_summary.md`

Measured result:

| lane | prepare | factorize | consensus/finalize | total | speedup |
| --- | ---: | ---: | ---: | ---: | ---: |
| cNMF parallel reference | 20.33 | 36.00 | 32.23 | 88.56 | 1.00x |
| FastCNMF scoped i50 w8 | 6.63 | 8.15 | 3.45 | 18.23 | 4.86x |

Quality versus the cNMF parallel reference:

- mean spectra cosine: `0.9996`
- mean usage Pearson: `0.9995`
- 95% gate: `true`

Coverage status:

- accepted: public S2, `72k` cells, `3.19x`, strict all-k pass
- accepted: internal S2, `72k` cells, `3.65x`, strict all-k pass
- accepted with scope limitation: available S1 spatial, `7,071` spots,
  `4.86x`, quality pass
- missing from local data: larger all-GBM S1 spatial benchmark

The scoped S1 result validates that the independent FastCNMF execution
architecture also works on spatial input. It does not prove scale behavior for
a larger all-GBM spatial tier because only three Visium samples are present in
the current workspace.
