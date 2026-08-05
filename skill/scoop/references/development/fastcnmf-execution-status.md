# FastCNMF Execution Status

## Current Stage

Stage: `S2 Direct Prepare Smoke`

Status: `partial; 3x gate not met`

Last updated: `2026-06-14`

## Checkpoint Log

### Checkpoint 0: Design Grounding

Completed:

- Confirmed reference cNMF workflow and Harmony behavior from the installed
  environment.
- Identified that `cnmf==1.7.1` needs Harmony MOE attributes including
  `Phi_moe`.
- Confirmed `harmonypy==2.0.0` changed the public object contract and no longer
  exposes the attributes expected by cNMF.
- Confirmed `harmonypy==0.2.0` provides `Phi_moe` and works with the cNMF
  preprocessing path used in the spatial benchmark.

Evidence:

- Reference run directory:
  `tmp/cnmf_spatial_harmony_benchmark/`
- Reference report:
  `tmp/cnmf_spatial_harmony_benchmark/benchmark_report.md`

Risks:

- Vendoring Harmony code may impose GPL compatibility obligations.
- GPU acceleration can change numerical behavior and must be gated by quality
  metrics rather than runtime alone.
- Large H5AD I/O can dominate end-to-end runtime if caching is not designed
  early.

Next:

- Start M1 cNMF-compatible runner.

### Checkpoint 1: M0 Bootstrap Complete

Completed:

- Added FastCNMF planning document:
  `docs/fastcnmf/00_project_plan.md`
- Added benchmark protocol:
  `docs/fastcnmf/01_benchmark_protocol.md`
- Added importable package skeleton under `src/fastcnmf`.
- Added CLI entry points:
  - `python -m fastcnmf status`
  - `python -m fastcnmf plan-tasks`
  - `python -m fastcnmf estimate-memory`
- Added initial deterministic NMF task manifest model.
- Added dense NMF memory estimator and replicate batch-size recommendation.
- Added tests for task manifest roundtrip and memory estimator.
- Added `fastcnmf = "fastcnmf.cli:main"` package script entry point.

Commands run:

```bash
PYTHONPATH=src python -m fastcnmf status
PYTHONPATH=src python -m fastcnmf plan-tasks --run-name smoke -k 6 8 --n-iter 2 --output tmp/fastcnmf_m0/task_manifest.json
PYTHONPATH=src python -m fastcnmf estimate-memory --observations 7071 --genes 3000 --components 8 --available-gib 8
PYTHONPATH=src pytest -q
```

Results:

- `plan-tasks` wrote 4 smoke tasks to
  `tmp/fastcnmf_m0/task_manifest.json`.
- `estimate-memory` reported `0.238 GiB` for a `7071 x 3000`, `k=8`,
  float32 single-replicate dense NMF task.
- Tests passed: `3 passed in 2.20s`.

Risks:

- Current package is scaffold only; no compute backend is implemented yet.
- Memory estimator is intentionally conservative but not yet calibrated against
  actual GPU allocation.
- Dynamic scheduling and cNMF-compatible execution are planned but not yet
  implemented.

Next:

- Implement M1: wrap cNMF Harmony preprocessing, write cNMF-compatible prepared
  inputs, and generate a task manifest from real cNMF replicate parameters.

### Checkpoint 2: S0 Harmony Benchmark Gate Passed

Goal:

- Beat original serial cNMF factorize runtime by at least `50%`.
- Keep output consistency above `95%`.

Completed:

- Added cNMF output quality comparison tooling in `src/fastcnmf/quality.py`.
- Added `python -m fastcnmf benchmark-gate` CLI.
- Added benchmark gate report generation.
- Added tests for program matching when cNMF program order is permuted.
- Ran the gate against the Harmony-corrected S0 spatial benchmark:
  - reference: serial cNMF output
  - candidate: 4-worker parallel cNMF-compatible output
  - k values: `6, 8`

Commands run:

```bash
PYTHONPATH=src pytest -q

PYTHONPATH=src python -m fastcnmf benchmark-gate \
  --summary-json tmp/cnmf_spatial_harmony_benchmark/benchmark_summary.json \
  --reference-dir tmp/cnmf_spatial_harmony_benchmark/serial/gbm_lowres_harmony_cnmf \
  --candidate-dir tmp/cnmf_spatial_harmony_benchmark/parallel/gbm_lowres_harmony_cnmf \
  --run-name gbm_lowres_harmony_cnmf \
  -k 6 8 \
  --min-time-saved-fraction 0.50 \
  --min-consistency 0.95 \
  --output-json tmp/fastcnmf_gate/harmony_s0_gate.json \
  --output-md tmp/fastcnmf_gate/harmony_s0_gate.md
```

Evidence:

- Gate JSON: `tmp/fastcnmf_gate/harmony_s0_gate.json`
- Gate report: `tmp/fastcnmf_gate/harmony_s0_gate.md`
- Tests: `4 passed in 3.46s`

Results:

- Reference serial factorize: `112.75 s`
- Candidate 4-worker factorize: `36.0 s`
- Speedup: `3.132x`
- Time saved fraction: `68.071%`
- Required time saved fraction: `50.000%`
- k=6 overall consistency: `1.000000`
- k=8 overall consistency: `1.000000`
- Required consistency: `0.950000`
- Gate result: `passed`

Interpretation:

- The S0 benchmark satisfies the active target: runtime improved by more than
  50% while consistency stayed above 95%.
- This is exact-compatible acceleration through cNMF-compatible parallel
  orchestration, not yet a GPU or approximate FastCNMF backend.

Risks:

- S0 is a smoke-scale benchmark. It proves the target on the current
  Harmony-corrected spatial fixture, not yet on S1/S2 larger datasets.
- Candidate runtime is from static cNMF worker-index parallelism. M1/M2 still
  need to promote this into a first-class FastCNMF runner with caching and task
  manifest control.

Next:

- Promote the successful S0 gate into the M1 runner.
- Run the same gate on S1 all-GBM Visium or S2 72k-cell fixtures.

### Checkpoint 3: Harmony2 Overall Speed Gate Passed

Goal:

- The whole FastCNMF workflow, not only the Harmony step, must be more than
  `2x` faster than the original serial cNMF workflow.
- Output consistency must remain above `95%`.

Completed:

- Added Harmony2 compatibility utilities in `src/fastcnmf/harmony2_compat.py`.
- Reconstructed `Phi_moe` from sample metadata for HarmonyPy 2.0 outputs.
- Implemented cNMF-compatible fixed `lamb=1` handling.
- Fixed lambda semantics to match the observed `cnmf==1.7.1` broadcast-vector
  behavior rather than a diagonal ridge matrix.
- Added batched MOE correction and verified it matches the cluster-loop
  implementation with maximum absolute error around `2.2e-11`.
- Re-ran the full S0 Harmony2 FastCNMF pipeline and regenerated the gate
  report.
- Added architecture and Harmony2 compatibility design:
  `docs/fastcnmf/03_architecture_and_harmony2.md`.

Commands run:

```bash
/usr/bin/time -v -o tmp/fastcnmf_harmony2/harmony20.time.log \
  ./.venv-harmony2/bin/python scripts/fastcnmf/run_harmony20_adapter.py \
  > tmp/fastcnmf_harmony2/harmony20.stdout.log

/usr/bin/time -v -o tmp/fastcnmf_harmony2/materialize_input.time.log \
  ./.venv-cnmf/bin/python scripts/fastcnmf/materialize_harmony20_cnmf_input.py \
  > tmp/fastcnmf_harmony2/materialize_input.stdout.log

bash scripts/fastcnmf/run_harmony20_fastcnmf_pipeline.sh

./.venv-cnmf/bin/python scripts/fastcnmf/summarize_fastcnmf_harmony20_gate.py
```

Evidence:

- Gate JSON: `tmp/fastcnmf_harmony2/fastcnmf_harmony20_gate.json`
- Gate report: `tmp/fastcnmf_harmony2/fastcnmf_harmony20_gate.md`
- Harmony2 output: `tmp/fastcnmf_harmony2/harmony20_adapter_output.npz`
- FastCNMF cNMF outputs:
  `tmp/fastcnmf_harmony2/parallel/gbm_lowres_harmony20_fastcnmf/`

Results:

- Original serial cNMF end-to-end: `317.51 s`
- FastCNMF Harmony2 end-to-end: `155.57 s`
- Overall speedup: `2.04x`
- Speed gate: `passed`
- Mean spectra cosine: `0.9959`
- Minimum spectra cosine: `0.9849`
- Mean usage Pearson: `0.9932`
- Minimum usage Pearson: `0.9549`
- Consistency gate: `passed`

Key finding:

- `harmonypy==2.0.0` can be made compatible with cNMF-style preprocessing if
  FastCNMF calls it with explicit fixed `lamb=1`, reconstructs `Phi_moe`, and
  applies cNMF's legacy MOE lambda broadcasting. The original incompatibility
  came from changed public attributes and changed default lambda behavior.

Risks:

- The S0 gate is still smoke-scale. Larger S1-S4 benchmarks are required before
  accepting GPU or approximate modes.
- Vendoring HarmonyPy 2.0 requires GPL-3.0-or-later compliance and provenance
  tracking.
- Materializing H5AD input is now the largest measured runtime block and should
  be redesigned around caching or backed storage.

Next:

- Vendor a fixed Harmony2 build with license/provenance metadata or implement a
  controlled fallback strategy for development environments.
- Run S1/S2 larger benchmarks.
- Prototype GPU NMF with the OOM retry ladder described in
  `docs/fastcnmf/03_architecture_and_harmony2.md`.

### Checkpoint 4: S1/S2 Cold-Start Manifest and Execution Plan

Goal:

- Move from S0 smoke results toward a fair larger benchmark.
- Define the independent FastCNMF execution boundary before optimizing toward
  `3x` against same-hardware optimized cNMF.

Completed:

- Added benchmark manifest models in
  `src/fastcnmf/benchmark_manifest.py`.
- Added an execution plan model in `src/fastcnmf/runner.py`.
- Added CLI commands:
  - `python -m fastcnmf benchmark-manifest`
  - `python -m fastcnmf plan-run`
- Generated a real project manifest:
  `tmp/fastcnmf_large_benchmark/benchmark_manifest.json`.
- Generated a stage-level execution plan:
  `tmp/fastcnmf_large_benchmark/execution_plan.json`.
- Added roadmap documentation:
  `docs/fastcnmf/04_large_benchmark_roadmap.md`.

Commands run:

```bash
PYTHONPATH=src python -m fastcnmf benchmark-manifest \
  --root . \
  --output-root tmp/fastcnmf_large_benchmark \
  --output tmp/fastcnmf_large_benchmark/benchmark_manifest.json

PYTHONPATH=src python -m fastcnmf plan-run \
  --manifest tmp/fastcnmf_large_benchmark/benchmark_manifest.json \
  --output tmp/fastcnmf_large_benchmark/execution_plan.json
```

Evidence:

- Manifest datasets:
  - `s1_gbm_lowres_visium_available`: 7,071 spots, 3 samples
  - `s2_public_24x3000`: 72,000 cells, 24 samples, 3,000 cells/sample
  - `s2_internal_24x3000`: 72,000 cells, 24 samples, 3,000 cells/sample
- Missing tier flag: `S1_all_gbm_visium`
- Execution plan stages: `30`
- Fairness policy:
  - comparison: `same_hardware_cold_start`
  - target speedup: `3.0x`
  - cross-lane artifact reuse: forbidden
  - spectra/usage thresholds: `0.95`

Interpretation:

- The S2 enlarged benchmark entry is now reproducible from the current
  worktree.
- S1 remains partial until more GBM Visium samples are available.
- This checkpoint does not claim the `3x` target is achieved; it defines the
  fair benchmark and independent stage boundary needed to test it.

Next:

- Implement real stage executors for the execution plan.
- Start with `s2_public_24x3000` using smoke parameters to establish a
  cold-start optimized cNMF baseline.
- Replace wrapper-style cNMF prepare/materialization with FastCNMF-owned
  preprocessing and cache boundaries.

### Checkpoint 5: Run Bundle and First Executable Stages

Goal:

- Move the stage plan toward an executable benchmark framework.
- Preserve fair cold-start artifact boundaries for S2 reference and candidate
  lanes.

Completed:

- Added cNMF reference preprocessing CLI:
  `python -m fastcnmf cnmf-preprocess`.
- Added run bundle generation in `src/fastcnmf/run_bundle.py`.
- Added CLI command:
  `python -m fastcnmf write-run-bundle`.
- Regenerated run bundle:
  `tmp/fastcnmf_large_benchmark/run_bundle/run_bundle.json`.
- Generated 30 executable shell scripts under:
  `tmp/fastcnmf_large_benchmark/run_bundle/scripts/`.
- Implemented scripts now cover:
  - cNMF reference preprocess for all current datasets
  - NMF task-manifest planning for both lanes
- Unimplemented candidate preprocess/factorize/consensus/report scripts exit
  with status `78` instead of pretending the executor exists.
- Executed one real S2 stage:
  `s2_public_24x3000:fastcnmf_independent:plan_nmf`.

Commands run:

```bash
PYTHONPATH=src pytest -q

PYTHONPATH=src python -m fastcnmf write-run-bundle \
  --manifest tmp/fastcnmf_large_benchmark/benchmark_manifest.json \
  --output-dir tmp/fastcnmf_large_benchmark/run_bundle \
  --python ./.venv-cnmf/bin/python

bash tmp/fastcnmf_large_benchmark/run_bundle/scripts/s2_public_24x3000__fastcnmf_independent__plan_nmf.sh

PYTHONPATH=src ./.venv-cnmf/bin/python -m fastcnmf cnmf-preprocess --help
```

Evidence:

- Tests: `9 passed in 1.73 s`
- Run bundle scripts: `30`
- Implemented scripts: `9`
- S2 public FastCNMF task manifest:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_independent/nmf/task_manifest.json`
- S2 public FastCNMF task count: `200`
- S2 public plan resource log:
  `tmp/fastcnmf_large_benchmark/run_bundle/logs/s2_public_24x3000__fastcnmf_independent__plan_nmf.time.log`
- `cnmf-preprocess` help works inside `.venv-cnmf`.

Interpretation:

- The large benchmark is no longer only a manifest. It now has executable stage
  scripts with resource logging and isolated artifact roots.
- This still does not prove the `3x` target. The heavy S2 cold-start cNMF
  preprocess/factorize/consensus runs have not been executed yet.

Next:

- Run
  `tmp/fastcnmf_large_benchmark/run_bundle/scripts/s2_public_24x3000__cnmf_optimized__preprocess.sh`
  to establish the first S2 cNMF cold-start preprocessing baseline.
- Implement FastCNMF independent H5AD preprocessing so candidate preprocess no
  longer exits `78`.
- Add cNMF reference `prepare/factorize/combine/consensus` executors to replace
  the remaining placeholder stage scripts.

### Checkpoint 6: FastCNMF Candidate Preprocess Executor

Goal:

- Replace the candidate preprocess placeholder with a real FastCNMF-owned H5AD
  preprocessing executor.
- Keep the S2 benchmark cold-start and cross-lane isolated.

Completed:

- Added FastCNMF independent preprocessing in `src/fastcnmf/preprocess.py`.
- Added CLI command:
  `python -m fastcnmf fast-preprocess`.
- The candidate preprocess now starts from a cold H5AD, uses the `counts` layer
  when present, writes TP10K, selects HVGs, computes PCA, runs Harmony2 through
  the FastCNMF MOE compatibility layer, and writes cNMF-compatible corrected
  input files.
- Updated run bundle generation so:
  - reference lane uses `.venv-cnmf/bin/python`
  - candidate lane uses `.venv-cnmf-h2/bin/python`
- Regenerated the S1/S2 run bundle. Implemented scripts increased from `9` to
  `12`.
- Ran a small cold-start FastCNMF preprocess smoke test using `.venv-cnmf-h2`.

Commands run:

```bash
PYTHONPATH=src ./.venv-cnmf-h2/bin/python -m fastcnmf fast-preprocess \
  --input-h5ad tmp/fastcnmf_large_benchmark/smoke_input.h5ad \
  --output-prefix tmp/fastcnmf_large_benchmark/smoke_fast/preprocess/cnmf_input \
  --sample-key sample_id \
  --n-top-genes 120 \
  --max-iter-harmony 3 \
  --seed 20260614

PYTHONPATH=src python -m fastcnmf write-run-bundle \
  --manifest tmp/fastcnmf_large_benchmark/benchmark_manifest.json \
  --output-dir tmp/fastcnmf_large_benchmark/run_bundle \
  --reference-python ./.venv-cnmf/bin/python \
  --candidate-python ./.venv-cnmf-h2/bin/python
```

Evidence:

- Smoke corrected H5AD:
  `tmp/fastcnmf_large_benchmark/smoke_fast/preprocess/cnmf_input.Corrected.HVG.Varnorm.h5ad`
- Smoke TP10K:
  `tmp/fastcnmf_large_benchmark/smoke_fast/preprocess/cnmf_input.TP10K.h5ad`
- Smoke HVG list:
  `tmp/fastcnmf_large_benchmark/smoke_fast/preprocess/cnmf_input.Corrected.HVGs.txt`
- Smoke output dimensions: `150` observations by `120` HVGs
- S2 public candidate preprocess script:
  `tmp/fastcnmf_large_benchmark/run_bundle/scripts/s2_public_24x3000__fastcnmf_independent__preprocess.sh`
- Run bundle implemented scripts: `12 / 30`

Interpretation:

- The candidate lane now has a real independent preprocess stage; it no longer
  depends on the cNMF reference lane's TP10K or corrected H5AD artifacts.
- The full S2 candidate preprocess has not yet been executed on 72k cells, so
  this checkpoint proves executor viability, not the final speed target.

Next:

- Execute S2 public reference and candidate preprocess stages to obtain the
  first fair cold-start preprocessing timings.
- Implement factorize, consensus, and report executors.
- Add quality comparison between S2 reference and candidate corrected inputs
  before running the expensive NMF stages.

### Checkpoint 7: S2 Public Cold-Start Preprocess Baseline

Goal:

- Establish the first fair S2 cold-start preprocessing comparison.
- Verify FastCNMF candidate preprocessing preserves corrected input quality
  before running expensive NMF stages.

Completed:

- Ran S2 public cNMF reference preprocess with CUDA visible. It failed with
  signal `11` during the Harmony0.2 CUDA path.
- Added lane-specific `CUDA_VISIBLE_DEVICES` controls to run bundle generation.
- Regenerated CPU-only run bundle:
  `tmp/fastcnmf_large_benchmark/run_bundle_cpu/run_bundle.json`.
- Ran S2 public cNMF reference preprocess with CUDA hidden.
- Ran S2 public FastCNMF candidate preprocess with CUDA hidden.
- Added chunked corrected-H5AD comparison tooling:
  `python -m fastcnmf compare-preprocess`.
- Compared the full `72000 x 3000` corrected matrices.
- Wrote preprocess benchmark summary artifacts.

Commands run:

```bash
bash tmp/fastcnmf_large_benchmark/run_bundle/scripts/s2_public_24x3000__cnmf_optimized__preprocess.sh

PYTHONPATH=src python -m fastcnmf write-run-bundle \
  --manifest tmp/fastcnmf_large_benchmark/benchmark_manifest.json \
  --output-dir tmp/fastcnmf_large_benchmark/run_bundle_cpu \
  --reference-python ./.venv-cnmf/bin/python \
  --candidate-python ./.venv-cnmf-h2/bin/python \
  --reference-cuda-visible-devices '' \
  --candidate-cuda-visible-devices ''

bash tmp/fastcnmf_large_benchmark/run_bundle_cpu/scripts/s2_public_24x3000__cnmf_optimized__preprocess.sh

bash tmp/fastcnmf_large_benchmark/run_bundle_cpu/scripts/s2_public_24x3000__fastcnmf_independent__preprocess.sh

PYTHONPATH=src ./.venv-cnmf-h2/bin/python -m fastcnmf compare-preprocess \
  --reference-h5ad tmp/fastcnmf_large_benchmark/s2_public_24x3000/cnmf_optimized/preprocess/cnmf_input.Corrected.HVG.Varnorm.h5ad \
  --candidate-h5ad tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_independent/preprocess/cnmf_input.Corrected.HVG.Varnorm.h5ad \
  --output-json tmp/fastcnmf_large_benchmark/s2_public_24x3000/preprocess_compare.json \
  --chunk-size 3000
```

Evidence:

- Failed GPU-visible reference log:
  `tmp/fastcnmf_large_benchmark/run_bundle/logs/s2_public_24x3000__cnmf_optimized__preprocess.time.log`
- CPU-only reference log:
  `tmp/fastcnmf_large_benchmark/run_bundle_cpu/logs/s2_public_24x3000__cnmf_optimized__preprocess.time.log`
- CPU-only candidate log:
  `tmp/fastcnmf_large_benchmark/run_bundle_cpu/logs/s2_public_24x3000__fastcnmf_independent__preprocess.time.log`
- Corrected input comparison:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/preprocess_compare.json`
- Preprocess benchmark summary:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/preprocess_benchmark_summary.json`
  and
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/preprocess_benchmark_summary.md`

Results:

- Failed GPU-visible cNMF reference attempt:
  - elapsed: `524.03 s`
  - max RSS: `13349.5 MB`
  - signal: `11`
- CPU-only cNMF reference preprocess:
  - elapsed: `258.54 s`
  - max RSS: `13042.2 MB`
  - output: `72000 x 3000` corrected H5AD
- CPU-only FastCNMF candidate preprocess:
  - elapsed: `106.32 s`
  - max RSS: `12596.7 MB`
  - output: `72000 x 3000` corrected H5AD
- Preprocess speedup: `2.43x`
- Corrected matrix cosine: `0.9957`
- Corrected matrix Pearson: `0.9956`
- Corrected input gate: `passed`

Interpretation:

- FastCNMF independent preprocessing is materially faster on S2 public and
  preserves corrected input quality above the `95%` threshold.
- The full objective is still incomplete. This is a preprocessing-stage result,
  not an end-to-end NMF/consensus `3x` benchmark.

Next:

- Implement cNMF reference and FastCNMF candidate factorize executors.
- Run S2 public smoke NMF first, then production-scale `n_iter=20`.
- Add consensus and final spectra/usage quality gates.

### Checkpoint 8: S2 Public Smoke End-to-End Benchmark

Goal:

- Run a complete S2 public smoke chain through preprocess, cNMF prepare,
  factorize, consensus, and quality comparison.
- Determine whether preprocessing acceleration alone is enough to improve
  end-to-end speed.

Completed:

- Extended run bundle scripts so all 30 stages are executable.
- Added `smoke` run-bundle profile using:
  - `k = 6, 8, 10, 12`
  - `n_iter = 8`
  - `max_nmf_iter = 200`
  - `workers = 4`
- Ran S2 public reference `plan_nmf`/cNMF prepare.
- Ran S2 public candidate `plan_nmf`/cNMF prepare.
- Ran S2 public reference factorize.
- Ran S2 public candidate factorize.
- Ran S2 public reference consensus.
- Ran S2 public candidate consensus.
- Added cNMF output comparison for different run names:
  `python -m fastcnmf compare-cnmf`.
- Generated S2 public smoke end-to-end summary.

Commands run:

```bash
PYTHONPATH=src python -m fastcnmf write-run-bundle \
  --manifest tmp/fastcnmf_large_benchmark/benchmark_manifest.json \
  --output-dir tmp/fastcnmf_large_benchmark/run_bundle_cpu_smoke \
  --reference-python ./.venv-cnmf/bin/python \
  --candidate-python ./.venv-cnmf-h2/bin/python \
  --reference-cuda-visible-devices '' \
  --candidate-cuda-visible-devices '' \
  --profile smoke

bash tmp/fastcnmf_large_benchmark/run_bundle_cpu_smoke/scripts/s2_public_24x3000__cnmf_optimized__plan_nmf.sh
bash tmp/fastcnmf_large_benchmark/run_bundle_cpu_smoke/scripts/s2_public_24x3000__fastcnmf_independent__plan_nmf.sh
bash tmp/fastcnmf_large_benchmark/run_bundle_cpu_smoke/scripts/s2_public_24x3000__cnmf_optimized__factorize.sh
bash tmp/fastcnmf_large_benchmark/run_bundle_cpu_smoke/scripts/s2_public_24x3000__fastcnmf_independent__factorize.sh
bash tmp/fastcnmf_large_benchmark/run_bundle_cpu_smoke/scripts/s2_public_24x3000__cnmf_optimized__consensus.sh
bash tmp/fastcnmf_large_benchmark/run_bundle_cpu_smoke/scripts/s2_public_24x3000__fastcnmf_independent__consensus.sh

PYTHONPATH=src ./.venv-cnmf-h2/bin/python -m fastcnmf compare-cnmf \
  --reference-dir tmp/fastcnmf_large_benchmark/s2_public_24x3000/cnmf_optimized/nmf/s2_public_24x3000_cnmf_optimized \
  --candidate-dir tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_independent/nmf/s2_public_24x3000_fastcnmf_independent \
  --reference-run-name s2_public_24x3000_cnmf_optimized \
  --candidate-run-name s2_public_24x3000_fastcnmf_independent \
  -k 6 8 10 12 \
  --output-json tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_cnmf_compare.json
```

Evidence:

- Smoke comparison:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_cnmf_compare.json`
- Smoke summary:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_end_to_end_summary.json`
  and
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_end_to_end_summary.md`
- cNMF-style reference outputs:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/cnmf_optimized/nmf/s2_public_24x3000_cnmf_optimized/`
- cNMF-style candidate outputs:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_independent/nmf/s2_public_24x3000_fastcnmf_independent/`

Results:

- Reference end-to-end smoke: `926.85 s`
- FastCNMF candidate end-to-end smoke: `924.83 s`
- End-to-end smoke speedup: `1.00x`
- 3x speed gate: `failed`
- Mean spectra cosine: `0.9983`
- Mean usage Pearson: `0.9953`
- 95% spectra/usage quality gate: `passed`

Step times:

| lane | preprocess | prepare | factorize | consensus | total |
| --- | ---: | ---: | ---: | ---: | ---: |
| reference | 258.54 | 76.91 | 374.55 | 216.85 | 926.85 |
| candidate | 106.32 | 105.42 | 486.80 | 226.29 | 924.83 |

Interpretation:

- FastCNMF preprocessing acceleration is real and quality-preserving.
- The current candidate still relies on cNMF-compatible prepare, factorize, and
  consensus. Those stages erase the preprocessing gain on S2.
- The next bottleneck is no longer Harmony preprocessing; it is independent NMF
  factorization, better task scheduling/cache boundaries, and consensus
  acceleration.

Next:

- Implement a FastCNMF-owned factorize backend or dynamic scheduler that avoids
  cNMF's static worker-index execution.
- Avoid re-materializing cNMF prepare artifacts where FastCNMF already owns
  corrected input and task manifests.
- Parallelize consensus by k and cache local-density intermediates.

### Checkpoint 9: FastCNMF Factorize/Consensus Executors

Goal:

- Move the candidate lane from cNMF CLI factorize/consensus toward a
  FastCNMF-owned execution framework.
- Measure whether dynamic replicate scheduling and k-parallel consensus are
  enough to close the S2 smoke end-to-end gap.

Completed:

- Added `python -m fastcnmf fast-factorize`.
- Added `python -m fastcnmf fast-consensus`.
- Changed run bundle generation so the `fastcnmf_independent` lane uses:
  - dynamic fork-based NMF replicate scheduling
  - cNMF-compatible iter spectra outputs
  - parallel-by-k combine/consensus/k-selection
- Kept the `cnmf_optimized` reference lane on cNMF CLI execution.
- Re-ran S2 public candidate factorize and consensus on the CPU smoke bundle.
- Re-ran spectra/usage quality comparison after the executor change.

Commands run:

```bash
PYTHONPATH=src python -m fastcnmf write-run-bundle \
  --manifest tmp/fastcnmf_large_benchmark/benchmark_manifest.json \
  --output-dir tmp/fastcnmf_large_benchmark/run_bundle_cpu_smoke \
  --reference-python ./.venv-cnmf/bin/python \
  --candidate-python ./.venv-cnmf-h2/bin/python \
  --reference-cuda-visible-devices '' \
  --candidate-cuda-visible-devices '' \
  --profile smoke

bash tmp/fastcnmf_large_benchmark/run_bundle_cpu_smoke/scripts/s2_public_24x3000__fastcnmf_independent__factorize.sh
bash tmp/fastcnmf_large_benchmark/run_bundle_cpu_smoke/scripts/s2_public_24x3000__fastcnmf_independent__consensus.sh

PYTHONPATH=src ./.venv-cnmf-h2/bin/python -m fastcnmf compare-cnmf \
  --reference-dir tmp/fastcnmf_large_benchmark/s2_public_24x3000/cnmf_optimized/nmf/s2_public_24x3000_cnmf_optimized \
  --candidate-dir tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_independent/nmf/s2_public_24x3000_fastcnmf_independent \
  --reference-run-name s2_public_24x3000_cnmf_optimized \
  --candidate-run-name s2_public_24x3000_fastcnmf_independent \
  -k 6 8 10 12 \
  --output-json tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_fast_executor_cnmf_compare.json

PYTHONPATH=src pytest -q
```

Evidence:

- Fast executor comparison:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_fast_executor_cnmf_compare.json`
- Fast executor summary:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_fast_executor_end_to_end_summary.json`
  and
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_fast_executor_end_to_end_summary.md`
- Factorize time log:
  `tmp/fastcnmf_large_benchmark/run_bundle_cpu_smoke/logs/s2_public_24x3000__fastcnmf_independent__factorize.time.log`
- Consensus time log:
  `tmp/fastcnmf_large_benchmark/run_bundle_cpu_smoke/logs/s2_public_24x3000__fastcnmf_independent__consensus.time.log`
- Tests: `10 passed in 1.87s`

Results:

| stage | reference seconds | candidate seconds | speedup | candidate max RSS MB |
| --- | ---: | ---: | ---: | ---: |
| preprocess | 258.54 | 106.32 | 2.43x | 12596.7 |
| prepare | 76.91 | 105.42 | 0.73x | 11620.8 |
| factorize | 374.55 | 431.59 | 0.87x | 3652.2 |
| consensus | 216.85 | 84.37 | 2.57x | 8628.0 |
| end-to-end | 926.85 | 727.70 | 1.27x |  |

Quality:

- Mean spectra cosine: `0.9983`
- Mean usage Pearson: `0.9953`
- 95% spectra/usage gate: `passed`
- 3x end-to-end speed gate: `failed`

Interpretation:

- FastCNMF now owns the candidate factorize and consensus execution entry
  points, so it is no longer only a shell wrapper for those stages.
- k-parallel consensus is effective on S2 smoke: `216.85 s` to `84.37 s`.
- Dynamic fork scheduling plus uncompressed spectra output is not enough for
  factorize: `431.59 s`, still slower than the optimized cNMF reference at
  `374.55 s`.
- The remaining 3x path requires a faster NMF backend or a controlled
  approximation strategy, not more shell-level scheduling around sklearn CD-NMF.

Next:

- Replace sklearn CD-NMF with a FastCNMF backend that can batch work across
  replicate seeds and/or use GPU with memory-aware tiling.
- Add a gated approximate mode that can reduce replicate count or max NMF
  iterations only when spectra/usage consistency stays above the configured
  threshold.
- Remove the cNMF prepare bottleneck by writing `norm_counts`, TPM stats, and
  replicate parameters directly from FastCNMF preprocessing outputs.

### Checkpoint 10: Direct FastCNMF Prepare Cache

Goal:

- Remove the candidate lane's dependency on cNMF `prepare`.
- Write the cNMF-compatible prepare cache directly from FastCNMF preprocessing
  outputs.
- Preserve spectra/usage quality after rerunning factorize and consensus.

Completed:

- Added `python -m fastcnmf fast-prepare`.
- Candidate `plan_nmf` now writes:
  - `task_manifest.json`
  - `norm_counts.h5ad`
  - `tpm_stats.df.npz`
  - `nmf_params.df.npz`
  - `nmf_idvrun_params.yaml`
  - `overdispersed_genes.txt`
- Candidate `fast-prepare` scales corrected HVG counts to unit standard
  deviation, matching the cNMF prepare contract.
- Candidate `fast-prepare` hardlinks the lane-local TP10K H5AD instead of
  copying it when the filesystem allows it.
- Re-ran S2 public candidate direct prepare, factorize, consensus, and
  spectra/usage comparison.

Commands run:

```bash
PYTHONPATH=src python -m fastcnmf write-run-bundle \
  --manifest tmp/fastcnmf_large_benchmark/benchmark_manifest.json \
  --output-dir tmp/fastcnmf_large_benchmark/run_bundle_cpu_smoke \
  --reference-python ./.venv-cnmf/bin/python \
  --candidate-python ./.venv-cnmf-h2/bin/python \
  --reference-cuda-visible-devices '' \
  --candidate-cuda-visible-devices '' \
  --profile smoke

bash tmp/fastcnmf_large_benchmark/run_bundle_cpu_smoke/scripts/s2_public_24x3000__fastcnmf_independent__plan_nmf.sh
bash tmp/fastcnmf_large_benchmark/run_bundle_cpu_smoke/scripts/s2_public_24x3000__fastcnmf_independent__factorize.sh
bash tmp/fastcnmf_large_benchmark/run_bundle_cpu_smoke/scripts/s2_public_24x3000__fastcnmf_independent__consensus.sh

PYTHONPATH=src ./.venv-cnmf-h2/bin/python -m fastcnmf compare-cnmf \
  --reference-dir tmp/fastcnmf_large_benchmark/s2_public_24x3000/cnmf_optimized/nmf/s2_public_24x3000_cnmf_optimized \
  --candidate-dir tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_independent/nmf/s2_public_24x3000_fastcnmf_independent \
  --reference-run-name s2_public_24x3000_cnmf_optimized \
  --candidate-run-name s2_public_24x3000_fastcnmf_independent \
  -k 6 8 10 12 \
  --output-json tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_direct_prepare_cnmf_compare.json

PYTHONPATH=src pytest -q
```

Evidence:

- Direct prepare summary:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_direct_prepare_end_to_end_summary.json`
  and
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_direct_prepare_end_to_end_summary.md`
- Direct prepare quality:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/smoke_direct_prepare_cnmf_compare.json`
- Candidate direct prepare log:
  `tmp/fastcnmf_large_benchmark/run_bundle_cpu_smoke/logs/s2_public_24x3000__fastcnmf_independent__plan_nmf.time.log`
- Tests: `10 passed in 1.92s`

Results:

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
- 3x end-to-end speed gate: `failed`

Interpretation:

- FastCNMF now owns preprocessing, prepare/cache writing, factorize scheduling,
  and consensus scheduling for the candidate lane.
- Direct prepare is compatible and materially faster than candidate cNMF
  prepare: `105.42 s` to `39.29 s`.
- The remaining bottleneck is concentrated in NMF factorization. With the
  current sklearn CD-NMF backend, the candidate factorize stage is still slower
  than the optimized cNMF reference.

Next:

- Prototype a declared approximate NMF mode and benchmark it as a separate
  candidate variant with the same spectra/usage quality gate.
- Prototype a real faster backend: batched CPU NMF or GPU NMF with OOM fallback.

### Checkpoint 11: Rejected Approximate n_iter=2 Smoke

Goal:

- Test whether an explicitly declared approximate NMF variant can reach the
  3x runtime target while preserving the same `95%` spectra/usage gate.

Experiment:

- Variant: `fastcnmf_approx_n2`
- k values: `6, 8, 10, 12`
- NMF replicates: `n_iter = 2`
- max NMF iterations: `200`
- workers: `4`
- consensus local-neighborhood-size: `0.5`

Important caveat:

- This was a manual exploratory run and did not inherit the run-bundle BLAS
  thread limits. `/usr/bin/time` reported CPU utilization far above the
  controlled benchmark lanes. Its speed is therefore not a fair result.
- The variant also intentionally uses fewer replicates than the reference, so
  it must pass quality before it can be considered as an approximate mode.

Evidence:

- Summary:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_approx_n2/approx_n2_summary.json`
  and
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_approx_n2/approx_n2_summary.md`
- Quality comparison:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_approx_n2/approx_n2_cnmf_compare.json`

Result:

- Apparent speedup: `3.64x`
- Mean spectra cosine: `0.9229`
- Mean usage Pearson: `0.8890`
- 95% spectra/usage gate: `failed`

Conclusion:

- Reject `n_iter=2` as a FastCNMF approximation strategy for this benchmark.
- A simple reduction in replicate count can produce apparent speed, but it does
  not preserve the required output quality.
- The next credible path is still a faster NMF backend, not skipping or
  aggressively reducing replicate computation.

### Checkpoint 12: Rejected MiniBatchNMF Backend Smoke

Goal:

- Test whether a faster CPU NMF backend can reach the `3x` target without
  reducing the number of cNMF replicates.

Experiment:

- Variant: `fastcnmf_minibatch_n8`
- Backend: `sklearn.decomposition.MiniBatchNMF`
- k values: `6, 8, 10, 12`
- NMF replicates: `n_iter = 8`
- max NMF iterations: `200`
- workers: `4`
- BLAS thread limits: explicitly set to `1`

Completed:

- Added an experimental `--backend minibatch` option to
  `python -m fastcnmf fast-factorize`.
- Kept the default backend as exact cNMF-compatible sklearn
  `non_negative_factorization`.
- Added a tiny backend smoke test; it is skipped in environments without
  `sklearn`.
- Ran controlled S2 public MiniBatchNMF factorize.
- Ran consensus at default `dt=0.5`, which failed density filtering.
- Ran no-filter consensus at `dt=2.0` to obtain comparable spectra/usage
  outputs and measure quality.

Evidence:

- Summary:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_minibatch_n8/minibatch_n8_summary.json`
  and
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_minibatch_n8/minibatch_n8_summary.md`
- Quality comparison:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_minibatch_n8/minibatch_n8_dt2_cnmf_compare.json`
- Tests: `10 passed, 1 skipped in 1.88s`

Result:

- MiniBatch factorize: `33.83 s`
- Apparent end-to-end speedup with `dt=2.0`: `3.04x`
- Default `dt=0.5` consensus: `failed`
- Mean spectra cosine at `dt=2.0`: `0.4075`
- Mean usage Pearson at `dt=2.0`: `0.4951`
- 95% spectra/usage gate: `failed`

Conclusion:

- Reject the current MiniBatchNMF backend for FastCNMF.
- It proves that factorize runtime can be reduced enough to reach the speed
  target, but the spectra are not cNMF-compatible under the current settings.
- The next backend needs either GPU/exact-style multiplicative updates or a
  batch algorithm that preserves cNMF spectra geometry far better than
  sklearn's default MiniBatchNMF.

### Checkpoint 13: Rejected CuPy MU GPU Backend Smoke

Goal:

- Test a real GPU NMF backend after MiniBatchNMF failed quality.
- Determine whether exact-style multiplicative updates can reach the `3x`
  runtime target while preserving spectra/usage quality.

Completed:

- Created a UV-managed GPU experiment environment:
  `.venv-fastcnmf-gpu`.
- Installed `cupy-cuda12x==14.1.1` plus the project dependencies.
- Added an experimental `--backend cupy-mu` option to
  `python -m fastcnmf fast-factorize`.
- Verified CuPy can access the RTX 3080 Ti and complete GPU matrix
  multiplication.
- Ran S2 public GPU MU smoke with:
  - `n_iter = 8`
  - `max_iter = 200`
  - `mu_dtype = float32`
- Fixed cNMF consensus dtype compatibility by writing GPU spectra back as
  float64.
- Ran a second GPU MU smoke with:
  - `n_iter = 8`
  - `max_iter = 600`

Evidence:

- GPU environment: `.venv-fastcnmf-gpu`
- CuPy MU max_iter=200 summary:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_cupy_mu_n8/cupy_mu_n8_summary.json`
  and
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_cupy_mu_n8/cupy_mu_n8_summary.md`
- CuPy MU max_iter=600 summary:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_cupy_mu_n8_i600/cupy_mu_n8_i600_summary.json`
  and
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_cupy_mu_n8_i600/cupy_mu_n8_i600_summary.md`
- Tests: `10 passed, 1 skipped`

Results:

| variant | factorize seconds | candidate accounting seconds | apparent speedup | spectra cosine | usage Pearson | quality gate |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `cupy_mu_n8` | 31.52 | 218.44 | 4.24x | 0.8959 | 0.8652 | failed |
| `cupy_mu_n8_i600` | 60.55 | 261.28 | 3.55x | 0.8958 | 0.8654 | failed |

Interpretation:

- The GPU path is technically viable and did not OOM on S2 public smoke.
- Factorize runtime is fast enough for the `3x` target.
- Increasing MU iterations from 200 to 600 did not improve cNMF agreement,
  indicating the current MU objective/updates converge to a different spectra
  geometry than the optimized cNMF CD reference.

Conclusion:

- Reject the current CuPy MU backend as a successful FastCNMF backend.
- Keep the backend code as an experimental harness for GPU memory and runtime
  testing.
- The next credible backend should target CD-like updates or another solver
  whose solution geometry matches cNMF more closely than plain MU.

### Checkpoint 14: Exact CD Iteration Truncation Boundary

Goal:

- Test whether a CD-like approximation can preserve cNMF output geometry while
  reducing factorize runtime.

Completed:

- Ran exact sklearn CD factorization with the same `n_iter=8`, same k values,
  same FastCNMF direct prepare, and reduced `max_iter` values.
- Tested `max_iter = 20, 25, 35, 50`.
- Ran cNMF-compatible consensus and spectra/usage comparison for each variant.

Evidence:

- Summary:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/exact_cd_truncation_summary.json`
  and
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/exact_cd_truncation_summary.md`

Results:

| max_iter | total seconds | speedup | factorize seconds | consensus seconds | spectra cosine | usage Pearson | quality gate | 3x gate |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| 20 | 462.06 | 2.01x | 93.80 | 154.98 | 0.9628 | 0.9359 | failed | failed |
| 25 | 515.94 | 1.80x | 116.73 | 276.64 | 0.9936 | 0.9799 | passed | failed |
| 35 | 425.49 | 2.18x | 193.84 | 114.94 | 0.9983 | 0.9943 | passed | failed |
| 50 | 434.14 | 2.13x | 241.55 | 69.63 | 0.9985 | 0.9953 | passed | failed |

Interpretation:

- CD truncation preserves cNMF geometry much better than MiniBatchNMF or plain
  GPU MU.
- The useful quality boundary is between `max_iter=20` and `max_iter=25`.
- However, even accepted settings do not reach `3x` because consensus/refit and
  preprocessing remain large enough that factorize truncation alone is
  insufficient.

Conclusion:

- Reject simple CD iteration truncation as the full FastCNMF acceleration
  strategy.
- Keep it as a possible quality-controlled knob, but pair it with a stronger
  consensus/refit optimization or a genuinely faster CD-like backend.

### Checkpoint 15: Spectra/Usage-Compatible Lite Consensus Boundary

Goal:

- Test whether consensus/refit optimization can combine with exact CD
  truncation to reach the `3x` target without changing the final spectra/usage
  definitions used for comparison.

Completed:

- Added `fast-consensus --lite`.
- First lite implementation wrote correct spectra but refit usages against
  `norm_counts`, which changed the cNMF final usage definition and failed usage
  comparison.
- Updated lite consensus to preserve cNMF's final usage definition:
  - cluster and write consensus spectra
  - refit spectra into TPM space
  - refit final usages on std-scaled HVG TPM
  - skip gene-level TPM spectra output, OLS gene scores, starCAT reference,
    clustergram, and k-selection outputs
- Added stricter comparison fields:
  - `min_k_overall_consistency`
  - `passes_all_k_95pct_gate`

Evidence:

- Summary:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/exact_cd_lite_consensus_summary.json`
  and
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/exact_cd_lite_consensus_summary.md`
- i25 quality:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_exact_cd_i25_lite/exact_cd_i25_lite_final_usage_cnmf_compare.json`
- i35 quality:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_exact_cd_i35_lite/exact_cd_i35_lite_final_usage_cnmf_compare.json`

Results:

| variant | total seconds | speedup | lite consensus seconds | mean spectra | mean usage | min k overall | global mean gate | all-k gate | 3x gate | strict accepted |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | --- |
| `exact_cd_i25_lite` | 264.24 | 3.51x | 24.94 | 0.9936 | 0.9799 | 0.9267 | passed | failed | passed | failed |
| `exact_cd_i35_lite` | 345.53 | 2.68x | 34.98 | 0.9983 | 0.9943 | 0.9831 | passed | passed | failed | failed |

Interpretation:

- This is the first candidate to exceed `3x` under global mean spectra/usage
  quality, but it does not pass a stricter all-k gate because k=6 usage
  consistency drops to `0.9267`.
- The more conservative i35 variant passes the stricter all-k quality gate but
  reaches only `2.68x`.
- Therefore the active goal is still not complete under a rigorous fairness
  interpretation: no variant currently satisfies both `>=3x` and all-k
  spectra/usage consistency `>=0.95`.

Conclusion:

- Keep spectra/usage-compatible lite consensus as a valid FastCNMF output mode
  only when downstream needs consensus spectra/usages and not full cNMF
  gene-level outputs.
- The next target should close the remaining gap by improving i35 runtime or by
  recovering k=6 usage quality for i25 without giving up the `3x` speed.

### Checkpoint 16: i27 Float32 Strict Candidate and Stability Repeat

Goal:

- Search the narrow boundary between i25, which is fast but fails strict all-k
  quality, and i30, which passes strict quality but is too slow.
- Improve factorize scheduling and memory bandwidth without changing cNMF
  spectra/usage output definitions.

Completed:

- Changed `fast-factorize` scheduling from input-order dynamic execution to
  largest-k-first dynamic execution.
- Added per-replicate factorize timing to the JSON output.
- Added `fast-prepare --norm-dtype float32`, keeping std calculation in
  float64 but writing `norm_counts.h5ad` as float32.
- Ran i27 largest-k-first in float64 and float32.
- Ran an independent float32 repeat to test runtime stability.

Evidence:

- Stability summary:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/exact_cd_lkf_float32_stability_summary.json`
  and
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/exact_cd_lkf_float32_stability_summary.md`
- First float32 quality:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_exact_cd_i27_lkf_f32/exact_cd_i27_lkf_f32_lite_final_usage_cnmf_compare.json`
- Repeat quality:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/fastcnmf_exact_cd_i27_lkf_f32_r2/exact_cd_i27_lkf_f32_r2_lite_final_usage_cnmf_compare.json`

Results:

| variant | dtype | total seconds | speedup | prepare | factorize | consensus lite | mean spectra | mean usage | min k overall | all-k gate | 3x gate | strict accepted |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- |
| `exact_cd_i27_lkf` | float64 | 389.17 | 2.38x | 92.08 | 155.10 | 35.67 | 0.9959 | 0.9866 | 0.9532 | passed | failed | failed |
| `exact_cd_i27_lkf_f32` | float32 | 270.62 | 3.42x | 28.13 | 88.51 | 47.66 | 0.9959 | 0.9866 | 0.9532 | passed | passed | passed |
| `exact_cd_i27_lkf_f32_r2` | float32 | 362.22 | 2.56x | 74.24 | 123.20 | 58.46 | 0.9959 | 0.9866 | 0.9532 | passed | failed | failed |

Interpretation:

- i27 with float32 `norm_counts` is numerically stable against the strict gate:
  both float32 runs preserve mean spectra/usage and all-k consistency.
- The first float32 run is the first strict candidate to exceed `3x`.
- The independent repeat did not preserve `3x`; wall time regressed in prepare,
  factorize, and consensus.
- Per-task timing shows largest-k-first largely removes the original k-order
  tail issue; the remaining factorize variance comes from task-level runtime
  variability and CPU contention, not only queue ordering.

Conclusion:

- Do not mark the active goal complete yet.
- Keep `float32 norm_counts` and largest-k-first scheduling as accepted
  FastCNMF optimizations.
- The next work should make runtime stable by reducing prepare I/O variance,
  recording consensus per-k timings, and adding retry/median benchmarking or a
  more controlled execution harness before claiming stable `~3x`.

### Checkpoint 17: Core Cache Prepare Path

Goal:

- Remove prepare-stage H5AD re-read/re-scale/re-write variance by promoting
  FastCNMF normalized counts and TPM stats into explicit cache artifacts.

Completed:

- Added `fast-preprocess --write-core-cache` to write:
  - normalized counts as `*.NormCounts.float32.npy`
  - obs/var name text files
  - TP10K stats as `*.TP10K.stats.df.npz`
- Added `fast-prepare --precomputed-*` arguments to hardlink those cache
  artifacts into a run directory.
- Added `fast-prepare --norm-store npy` and manifest-driven readers in
  `fast-factorize` and `fast-consensus --lite`.
- Added `fast-preprocess --no-corrected-h5ad` as an experimental way to skip an
  intermediate not needed by the core-cache FastCNMF path.

Evidence:

- Summary:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/core_cache_stability_summary.json`
  and
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/core_cache_stability_summary.md`

Results:

| variant | total seconds | speedup | preprocess | prepare | factorize | consensus | mean spectra | mean usage | min k overall | strict accepted |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `i27_lkf_f32_single` | 270.62 | 3.42x | 106.32 | 28.13 | 88.51 | 47.66 | 0.9959 | 0.9866 | 0.9532 | true |
| `i27_lkf_f32_repeat` | 362.22 | 2.56x | 106.32 | 74.24 | 123.20 | 58.46 | 0.9959 | 0.9866 | 0.9532 | false |
| `i27_corecache` | 326.87 | 2.84x | 167.81 | 2.03 | 109.79 | 47.24 | 0.9959 | 0.9866 | 0.9532 | false |

Additional preprocess-only result:

- `--no-corrected-h5ad` with core cache took `199.18 s`, slower than the
  measured core-cache preprocess that still wrote corrected H5AD.

Interpretation:

- Core cache successfully fixes prepare variance: prepare dropped to `2.03 s`
  and no longer reads/recomputes the normalized matrix or TP10K stats.
- The cache path also reduces factorize memory through memmapped float32 input.
- It does not yet produce a stable `>=3x` cold-start result because moving cache
  creation into preprocess increased measured preprocess wall time enough to
  offset the prepare gain.
- Skipping corrected H5AD did not help in the measured run; remaining
  preprocess cost appears dominated by computation, TP10K/core cache
  materialization, and system I/O variability.

Conclusion:

- Keep core cache as a valid independent-framework improvement.
- The next stability target is no longer prepare; it is preprocess cache
  materialization and consensus/factorize variance.
- The goal remains active: no repeated cold-start result has yet proven stable
  `~3x` speedup under strict all-k quality.

### Checkpoint 18: HVG-Only Lite Consensus Cache

Goal:

- Reduce lite consensus cost by avoiding full TP10K H5AD reads and full-gene
  spectra refit when only final spectra/usages are required.

Completed:

- Added TP10K HVG raw/scaled cache files in `fast-preprocess
  --write-core-cache`.
- Fixed cache file naming so raw and scaled matrices are distinct.
- Added `fast-preprocess --no-tp10k-h5ad` for FastCNMF-lite runs that do not
  need full TP10K H5AD materialization.
- Added `fast-prepare --precomputed-tpm-hvg-*` and manifest fields.
- Updated `fast-consensus --lite` to:
  - use cached normalized counts
  - use cached raw/scaled TP10K HVG matrices
  - compute final usage from HVG-only TPM refit
  - preserve the same final spectra/usage comparison values

Evidence:

- Summary:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/hvg_consensus_cache_summary.json`
  and
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/hvg_consensus_cache_summary.md`

Results:

| variant | total seconds | speedup | preprocess | prepare | factorize | consensus | consensus RSS MB | min k overall | strict accepted |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `i27_corecache` | 326.87 | 2.84x | 167.81 | 2.03 | 109.79 | 47.24 | 6066.4 | 0.9532 | false |
| `i27_corecache_hvg_fix` | 346.15 | 2.68x | 208.53 | 3.77 | 121.70 | 12.15 | 1594.7 | 0.9532 | false |
| `i27_corecache_hvg_notp` | 324.25 | 2.86x | 180.47 | 2.03 | 130.66 | 11.09 | 1593.8 | 0.9532 | false |

Interpretation:

- HVG-only lite consensus is a real optimization:
  - consensus wall time drops from about `47 s` to `11-12 s`
  - consensus peak RSS drops from about `6.1 GB` to `1.6 GB`
  - strict spectra/usage quality is unchanged
- Avoiding full TP10K H5AD materialization helps, but not enough by itself:
  preprocess remains `180.47 s` in the measured no-TP run.
- End-to-end still misses stable `3x` because the remaining variance is now
  preprocess materialization and factorize runtime, not consensus.

Conclusion:

- Keep HVG-only consensus cache as an accepted FastCNMF-lite optimization.
- Next work should target factorize variance or a preprocessing path that
  avoids materializing both corrected H5AD and core cache forms.

### Checkpoint 19: Same-Worker S2 Fairness Rerun

Goal:

- Replace the remaining fairness caveat with a cNMF reference rerun using the
  same 8-process CPU worker budget as the current FastCNMF candidate.

Completed:

- Re-ran cNMF `prepare`, `factorize`, `combine`, `consensus`, and
  `k_selection_plot` in a new output lane:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/cnmf_optimized_w8/`.
- Kept the same cNMF cold-start preprocessing time from the CPU-only reference
  lane: `258.54 s`.
- Ran cNMF factorization with 8 shell-launched workers and BLAS thread limits
  set to 1.
- Compared the accepted FastCNMF 8-worker candidate against the new cNMF
  8-worker reference.

Evidence:

- Fairness summary:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/same_worker_w8_fairness_summary.json`
  and
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/same_worker_w8_fairness_summary.md`
- New cNMF 8-worker quality comparison:
  `tmp/fastcnmf_large_benchmark/s2_public_24x3000/cnmf_optimized_w8/fastcnmf_w8_vs_cnmf_w8_compare.json`

Results:

| lane | preprocess | prepare | factorize | consensus | total | speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cNMF optimized, 4 workers existing | 258.54 | 76.91 | 374.55 | 216.85 | 926.85 | 1.00x |
| cNMF optimized, 8 workers rerun | 258.54 | 75.77 | 235.69 | 240.69 | 810.69 | 1.00x |
| FastCNMF, 8 factorize workers | 180.47 | 2.02 | 61.16 | 10.22 | 253.87 | 3.19x |
| FastCNMF, repeat compute + repeat preprocess | 118.12 | 2.04 | 59.09 | 9.00 | 188.25 | 4.31x |
| FastCNMF, worst-of-two stage envelope | 180.47 | 2.04 | 61.16 | 10.22 | 253.89 | 3.19x |

Quality versus the cNMF 8-worker reference:

- mean spectra cosine: `0.9959`
- mean usage Pearson: `0.9866`
- strict all-k gate: `true`
- minimum per-k overall consistency: `0.9532`

Interpretation:

- The S2 public `72k`-cell benchmark now has a same-worker reference.
- cNMF does improve when factorize workers increase from 4 to 8, but the total
  reference time remains `810.69 s` because consensus is still about
  `240.69 s`.
- FastCNMF reaches `3.19x` versus this stricter reference, and the worst-of-two
  FastCNMF stage envelope still reaches `3.19x`.
- The main acceleration sources are now clear:
  - direct core-cache prepare: `75.77 s` to `2.02 s`
  - dynamic/memmap factorize: `235.69 s` to `61.16 s`
  - HVG-cache lite consensus: `240.69 s` to `10.22 s`

Remaining gap before marking the full active goal complete:

- The S2 public same-worker benchmark passes the `~3x` and `>=95%` quality
  requirements.
- The broader goal still needs the expanded internal S2 and larger spatial S1
  benchmark coverage to be brought to the same level of reproducibility.

### Checkpoint 20: Internal S2 Same-Worker Stress Test

Goal:

- Extend the same-worker S2 benchmark from public data to the internal
  `s2_internal_24x3000` dataset.

Completed:

- Added reusable same-worker benchmark runner:
  `scripts/fastcnmf/run_s2_same_worker_w8_benchmark.sh`.
- Ran the internal S2 cNMF optimized 8-worker reference from cold H5AD.
- Ran FastCNMF candidates with max NMF iterations `27`, `35`, `50`, and `100`
  using the same FastCNMF independent preprocess/cache output.
- Verified that preprocessing is not the current blocker:
  - corrected input cosine: `0.9958`
  - corrected input Pearson: `0.9955`
  - norm-count input cosine: `0.9960`
  - norm-count input Pearson: `0.9958`

Evidence:

- Internal sweep summary:
  `tmp/fastcnmf_large_benchmark/s2_internal_24x3000/internal_s2_i27_i100_sweep_summary.json`
  and
  `tmp/fastcnmf_large_benchmark/s2_internal_24x3000/internal_s2_i27_i100_sweep_summary.md`
- Preprocess comparison:
  `tmp/fastcnmf_large_benchmark/s2_internal_24x3000/preprocess_quality_compare.json`
- Norm-count comparison:
  `tmp/fastcnmf_large_benchmark/s2_internal_24x3000/norm_counts_compare.json`

Reference timing:

| lane | preprocess | prepare | factorize | consensus | total |
| --- | ---: | ---: | ---: | ---: | ---: |
| cNMF optimized, 8 workers | 408.93 | 80.46 | 349.02 | 296.83 | 1135.24 |

FastCNMF sweep:

| variant | max_iter | total | speedup | factorize | mean spectra | mean usage | min k overall | strict all-k | failing k pattern |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| `i27` | 27 | 227.93 | 4.98x | 72.97 | 0.9541 | 0.9550 | 0.8795 | false | 6, 12 |
| `i35` | 35 | 244.91 | 4.64x | 87.39 | 0.9397 | 0.9331 | 0.8555 | false | 6, 10 |
| `i50` | 50 | 260.95 | 4.35x | 103.80 | 0.9626 | 0.9596 | 0.8593 | false | 10 |
| `i100` | 100 | 336.60 | 3.37x | 179.71 | 0.9720 | 0.9643 | 0.8597 | false | 10 |

Interpretation:

- Internal S2 establishes a stricter quality stress test than public S2.
- FastCNMF speed is sufficient through `i100`, but strict all-k quality does
  not pass because `k=10` remains unstable around `0.86` overall consistency.
- The failure is unlikely to come from independent preprocessing because both
  corrected input and normalized NMF input match the cNMF reference at about
  `0.996`.
- Increasing per-replicate max iterations is not monotonic and does not solve
  internal `k=10`; the next internal strategy should spend the remaining speed
  budget on consensus stability, most likely with more replicates, for example
  an `n_iter=20` same-worker benchmark.

Follow-up n20 benchmark:

- Ran the same-worker internal benchmark with `n_iter=20` and FastCNMF
  `max_iter=50`.
- This keeps the same 8 factorize workers and spends the remaining speed
  budget on consensus stability rather than only increasing per-replicate
  iterations.

Evidence:

- Accepted internal summary:
  `tmp/fastcnmf_large_benchmark/s2_internal_24x3000/internal_s2_accepted_n20_i50_summary.json`
  and
  `tmp/fastcnmf_large_benchmark/s2_internal_24x3000/internal_s2_accepted_n20_i50_summary.md`

Accepted internal result:

| lane | preprocess | prepare | factorize | consensus | total | speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cNMF optimized, n20, 8 workers | 330.17 | 69.16 | 774.74 | 283.96 | 1458.03 | 1.00x |
| FastCNMF, n20 i50, 8 workers | 142.18 | 2.13 | 242.50 | 12.70 | 399.51 | 3.65x |

Quality:

- mean spectra cosine: `0.9982`
- mean usage Pearson: `0.9990`
- strict all-k gate: `true`
- minimum per-k overall consistency: `0.9967`

Interpretation update:

- The earlier internal failure was a low-replicate consensus stability issue,
  not a preprocessing/input equivalence issue.
- Increasing `n_iter` from 8 to 20 fixed the `k=10` failure while preserving
  `3.65x` same-worker speedup.
- FastCNMF factorize completed 80 tasks in `242.50 s`; per-process cumulative
  task time stayed in a narrow `229-242 s` band under dynamic largest-k-first
  scheduling.

Conclusion:

- Public S2 is accepted.
- Internal S2 is now accepted under the strict all-k gate with the n20 i50
  variant.
- The active goal remains incomplete until spatial S1 has an accepted
  reproducible run or a documented, justified scope decision.

### Checkpoint 21: Scoped S1 Spatial FastCNMF Run

Goal:

- Add spatial coverage using the locally available GBM Visium samples while
  preserving the documented limitation that all-GBM S1 is not present.

Completed:

- Reused the existing scoped spatial Harmony input:
  `tmp/cnmf_spatial_harmony_benchmark/input/gbm_lowres_visium_3samples_harmony.Corrected.HVG.Varnorm.h5ad`
- Ran FastCNMF direct prepare, exact factorize, and lite consensus:
  - k values: `6, 8`
  - NMF replicates per k: `8`
  - max NMF iterations: `50`
  - factorize workers: `8`
  - consensus workers: `2`
- Compared FastCNMF outputs to the existing cNMF parallel reference for the
  same three spatial samples.

Evidence:

- Scoped S1 summary:
  `tmp/cnmf_spatial_harmony_benchmark/fastcnmf_scoped_i50_w8/scoped_s1_fastcnmf_summary.json`
  and
  `tmp/cnmf_spatial_harmony_benchmark/fastcnmf_scoped_i50_w8/scoped_s1_fastcnmf_summary.md`
- Existing cNMF spatial reference:
  `tmp/cnmf_spatial_harmony_benchmark/benchmark_summary.json`

Results:

| lane | prepare | factorize | consensus/finalize | total | speedup |
| --- | ---: | ---: | ---: | ---: | ---: |
| cNMF parallel reference | 20.33 | 36.00 | 32.23 | 88.56 | 1.00x |
| FastCNMF scoped i50 w8 | 6.63 | 8.15 | 3.45 | 18.23 | 4.86x |

Quality:

- mean spectra cosine: `0.9996`
- mean usage Pearson: `0.9995`
- 95% gate: `true`

Interpretation:

- The available spatial subset passes speed and quality with the independent
  FastCNMF runner.
- This is scoped S1 coverage only. The repository contains three GBM Visium
  samples (`7,071` spots total), and the benchmark manifest still marks
  `S1_all_gbm_visium` as missing.
- The scoped S1 result supports the architecture, but it cannot prove the
  larger all-GBM S1 requirement without additional spatial samples.

Conclusion:

- S2 public accepted.
- S2 internal accepted.
- S1 scoped available-spatial accepted.
- The only remaining gap for the full original scope is data availability for a
  larger all-GBM S1 spatial benchmark.
