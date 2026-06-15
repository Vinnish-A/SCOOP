# FastCNMF Project Plan

## Objective

FastCNMF is a cNMF-compatible acceleration layer for large single-cell and
spatial transcriptomics program discovery. It keeps the cNMF/Harmony workflow
as the reference behavior while adding deterministic orchestration, caching,
parallel scheduling, GPU-capable NMF backends, and rigorous benchmark reporting.

The initial target is compatibility with the tested reference stack:

- `cnmf==1.7.1`
- `harmonypy==0.2.0` compatible Harmony MOE ridge correction
- AnnData/H5AD inputs with raw counts, `sample_id`, and optional spatial
  coordinates

## Non-Goals

- Do not replace cNMF consensus semantics in the compatibility mode.
- Do not require GPU for the base workflow.
- Do not silently use approximate algorithms when the user requested exact
  compatibility.
- Do not vendor Harmony code without carrying its license and provenance.

## Architecture

```text
fastcnmf/
  cli.py
  config.py
  resources.py

  io/
    h5ad.py
    zarr_store.py
    cnmf_export.py

  preprocess/
    harmony_compat.py
    hvg.py
    normalize.py

  nmf/
    backend.py
    sklearn_cpu.py
    torch_gpu.py
    minibatch.py

  scheduler/
    local.py
    task_manifest.py

  consensus/
    merge.py
    metrics.py

  benchmark/
    protocol.py
    report.py
```

## Execution Modes

### `exact-compatible`

Uses cNMF-compatible preprocessing, NMF parameters, task outputs, and consensus
format. Acceleration comes from orchestration, caching, resource control, and
parallel worker management.

Acceptance targets:

- spectra cosine similarity versus cNMF baseline: `>= 0.995`
- usage correlation versus cNMF baseline: `>= 0.99`
- output files readable by existing cNMF downstream tools

### `fast-compatible`

Uses float32 and optionally GPU NMF while keeping the same input preprocessing
and consensus contract.

Acceptance targets:

- spectra cosine similarity versus cNMF baseline: `>= 0.98`
- usage correlation versus cNMF baseline: `>= 0.95`
- reconstruction error delta: `<= 2%`

### `large-approx`

Uses sketching, mini-batch NMF, adaptive replicate stopping, and full-data
usage refit to support data that cNMF cannot practically complete.

Acceptance targets:

- spectra cosine similarity on validation subset: `>= 0.95`
- reconstruction error delta on validation subset: `<= 5%`
- no stronger `sample_id` leakage than cNMF/Harmony baseline

## Acceleration Points

| cNMF Bottleneck | FastCNMF Intervention |
| --- | --- |
| Harmony/correction dense materialization | block-wise MOE correction and optional Zarr backing |
| static worker-index scheduling | task manifest with dynamic local scheduler |
| repeated prepare and H5AD copying | prepared input cache and checksum validation |
| CPU-only replicate factorization | pluggable CPU/GPU NMF backends |
| consensus per-k serial execution | per-k parallel consensus tasks |
| failure requires manual cleanup | checkpointed task states and retry policy |

## GPU Memory Strategy

Before a GPU NMF task starts, FastCNMF estimates memory from:

- matrix chunk shape and dtype
- `W` and `H` factor shapes
- number of replicate tasks batched on the same device
- temporary buffers required by the selected NMF algorithm
- a configurable safety margin

If the estimate exceeds available VRAM, or a CUDA OOM is observed, the runtime
must retry with:

1. lower replicate batch size
2. smaller observation chunks
3. smaller gene chunks
4. float32 or mixed precision, if allowed by mode
5. another GPU device, if available
6. CPU backend fallback for the failed task

## Harmony Compatibility Plan

FastCNMF will provide a Harmony compatibility layer that exposes the quantities
cNMF requires for MOE ridge correction:

- `Z_corr`
- `R`
- `K`
- `Phi_moe`
- `lamb`

Short-term, FastCNMF pins and validates `harmonypy==0.2.0` behavior. Mid-term,
it can adapt `harmonypy==2.x` by reconstructing `Phi_moe` from metadata and
exposing final lambda values from the C++ backend or an equivalent
`apply_moe_correction` API.

## Milestones

### M0: Project Bootstrap

- Write project plan and benchmark protocol.
- Add importable `fastcnmf` package skeleton.
- Add CLI stubs for `preprocess`, `factorize`, `benchmark`, and `status`.
- Add status checkpoint document.

### M1: cNMF-Compatible Runner

- Wrap cNMF Harmony preprocessing with pinned compatibility.
- Generate cNMF-compatible prepared files.
- Create task manifest for NMF replicates.
- Run local static and dynamic worker benchmarks.

### M2: Benchmark Harness

- Run cNMF baseline and FastCNMF exact-compatible benchmark on S0/S1 data.
- Record wall time, CPU, RSS, I/O, and output quality metrics.
- Add automated report generation.

### M3: Dynamic Scheduler

- Replace static worker-index orchestration with task queue scheduling.
- Support retry and `skip_completed`.
- Compare dynamic versus cNMF static worker distribution.

### M4: GPU Backend Prototype

- Implement Torch/CuPy NMF backend behind a common interface.
- Add VRAM estimator and OOM retry ladder.
- Benchmark GPU on S1/S2 and record failure behavior.

### M5: Large-Scale Approximation

- Add sketching and full-data usage refit.
- Benchmark S3/S4 datasets.
- Quantify quality/runtime tradeoffs.

## Stage Gate Criteria

Every milestone must update `docs/fastcnmf/02_execution_status.md` with:

- completed work
- commands run
- outputs produced
- measured runtime/resource evidence
- risks and blockers
- next milestone decision

