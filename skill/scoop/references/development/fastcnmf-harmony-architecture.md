# FastCNMF Architecture and Harmony2 Compatibility

## Goal

FastCNMF is a cNMF-compatible acceleration layer. The first hard gate is not a
single-step microbenchmark: the whole workflow must be at least `2x` faster
than the original serial cNMF workflow while preserving `>=95%` output
consistency.

The current S0 spatial gate passes:

- original serial cNMF end-to-end: `317.51 s`
- FastCNMF Harmony2 end-to-end: `155.57 s`
- overall speedup: `2.04x`
- mean spectra cosine: `0.9959`
- mean usage Pearson: `0.9932`
- minimum usage Pearson: `0.9549`

Evidence:

- JSON: `tmp/fastcnmf_harmony2/fastcnmf_harmony20_gate.json`
- Markdown: `tmp/fastcnmf_harmony2/fastcnmf_harmony20_gate.md`

## Runtime Architecture

```text
input AnnData
  |
  v
FastCNMF preprocess cache
  - HVG selection
  - TP10K materialization
  - scaled HVG matrix
  - PCA
  |
  v
Harmony compatibility layer
  - vendored/pinned Harmony2 runner
  - reconstructed Phi_moe
  - cNMF-compatible lambda vector handling
  - batched MOE gene-space correction
  |
  v
cNMF-compatible prepared files
  - Corrected.HVG.Varnorm.h5ad
  - TP10K.h5ad
  - Corrected.HVGs.txt
  |
  v
NMF task scheduler
  - one task per k/replicate
  - static cNMF worker mode for exact compatibility
  - dynamic queue mode for better load balancing
  |
  v
consensus and export
  - cNMF-compatible spectra/usages
  - quality gate report
  - resource logs
```

## Acceleration Points

| Area | Current FastCNMF implementation | Next acceleration target |
| --- | --- | --- |
| Harmony execution | HarmonyPy 2.0 C++ backend | vendor fixed Harmony2 build and expose stable adapter |
| MOE correction | batched BLAS implementation, equivalent to cluster loop | block-wise correction for larger matrices and lower peak RAM |
| preprocessing | core arrays cached separately from cNMF files | avoid repeated H5AD writes with Zarr/backed arrays |
| factorization | cNMF worker-index parallelism, 4 workers in S0 | dynamic task queue to reduce idle workers |
| consensus | cNMF-compatible consensus | per-k parallel consensus and cached local density |
| GPU | memory estimator only | Torch/CuPy backend with OOM retry ladder |

The measured S0 bottlenecks after the Harmony2 fix are:

- materialize input: `49.93 s`
- factorize parallel: `37 s`
- core prepare: `17.85 s`
- consensus: `14.26 s`
- Harmony2 adapter process wall time: `11.36 s`

The largest remaining speed opportunity is reducing H5AD materialization and
moving from static worker-index scheduling to a dynamic task queue.

## Approximation Levels

FastCNMF should make approximation explicit. It must not silently trade quality
for speed.

| Mode | Contract | Minimum consistency | Intended dataset size |
| --- | --- | ---: | --- |
| `exact-compatible` | same cNMF output contract, cNMF-compatible Harmony MOE | spectra `>=0.995`, usage `>=0.99` | smoke and production baselines |
| `fast-compatible` | float32 and optional GPU, same preprocessing and consensus contract | spectra `>=0.98`, usage `>=0.95` | 50k-500k cells/spots |
| `large-approx` | sketch/refit, minibatch NMF, adaptive replicate stopping | validation spectra `>=0.95`, reconstruction delta `<=5%` | 500k-1M+ cells/spots |

Approximation should buy larger benchmarks, not just better timing on S0. A
reasonable large benchmark ladder is:

- S1: all available GBM Visium samples
- S2: public plus internal single-cell fixtures
- S3: 200k-500k sampled cells/spots
- S4: 1M cells/spots using sketch/refit

## GPU and OOM Policy

Before launching a GPU task, FastCNMF estimates memory for:

- dense input block
- `W` and `H`
- per-replicate temporary buffers
- concurrent replicate batch size
- dtype
- safety margin

If the estimate exceeds the available VRAM, FastCNMF should not launch that
batch. If CUDA still raises OOM, the task-level retry ladder is:

1. empty CUDA cache and retry once with the same task
2. reduce concurrent replicate batch size
3. reduce observation chunk size
4. reduce gene chunk size
5. switch to float32 or mixed precision if the selected mode permits it
6. move to another GPU if available
7. fall back to CPU for that task and mark the fallback in the report

The OOM benchmark must record estimated VRAM, actual peak allocated/reserved
VRAM, task size, retry count, and final backend. A GPU backend is not accepted
until it shows that OOM does not corrupt completed replicate outputs.

## HarmonyPy 0.2 versus 2.0

`harmonypy==0.2.0` works with `cnmf==1.7.1` because cNMF expects the Harmony
result object to expose:

- `Phi_moe`
- `lamb`
- `R`
- `K`
- `Z_corr`

It also matches cNMF's fixed-lambda MOE behavior. In the tested cNMF path,
`harmonypy 0.2` returns `lamb` as a vector like `[0, 1, 1, 1]`, and
`cnmf.preprocess.moe_correct_ridge` adds that vector directly to the Gram
matrix. This NumPy broadcast behavior is part of the observed compatibility
contract, even though it is not the diagonal ridge form one might expect.

`harmonypy==2.0.0` initially failed for two reasons:

- it no longer exposes `Phi_moe` and `lamb` on the public Harmony object
- its `lamb=None` behavior is dynamic lambda estimation, while the cNMF
  compatibility path needs fixed `lamb=1`

FastCNMF can use HarmonyPy 2.0 and keep cNMF compatibility by:

- calling Harmony2 with explicit `lamb=1`
- reconstructing `Phi_moe` from `obs[batch_key]`
- building the cNMF-compatible lambda vector
- applying MOE correction with the same broadcast-lambda semantics as cNMF
- using a batched implementation that matches the cluster loop numerically

The current batched MOE implementation matches the cluster-loop implementation
with maximum absolute error around `2.2e-11` on the S0 input.

## Vendoring Plan

HarmonyPy 2.0 is GPL-3.0-or-later. If FastCNMF vendors it, the package must
carry the Harmony license, copyright notices, and source provenance. The
recommended package layout is:

```text
src/fastcnmf/vendor/
  harmonypy2/
    harmonypy/
    harmonypy-2.0.0.dist-info/
    LICENSE
    PROVENANCE.md
```

FastCNMF code should import through a wrapper, not directly from external
`harmonypy`:

```python
from fastcnmf.harmony2_compat import harmony2_moe_correct
```

That wrapper becomes the stable contract. It can prefer the vendored fixed
Harmony2 build and optionally fall back to an installed `harmonypy==2.0.0` in
developer mode. Production benchmark reports must record the exact source:
vendored build hash or external wheel version.

## Benchmark Requirements

Every accepted FastCNMF acceleration must compare against upstream cNMF with:

- same input cells/spots and genes
- same `sample_id` batch-removal target
- same `k`, `n_iter`, seed, and max NMF iterations
- end-to-end wall clock, not just Harmony or NMF microbenchmarks
- CPU percent, max RSS, file I/O, and GPU metrics when applicable
- matched spectra cosine and usage Pearson by Hungarian assignment
- reconstruction error delta for approximate/GPU modes

The current S0 gate is a smoke benchmark. It proves the Harmony2 adapter and
overall `>2x` speedup on the selected GBM low-resolution spatial samples. It
does not replace S1-S4 larger benchmarks.
