# FastDE Abundance scSurvival Benchmark

FastDE abundance now uses a scSurvival-style multiple-instance learning backend
as the only active implementation. Earlier linear/CLR aggregation experiments
were removed because they did not preserve the core scSurvival idea that each
sample is a bag of cell instances.

## Reference Comparison Scope

The reference benchmark compares:

- original public `cliffren/scSurvival` at commit
  `a76de0a00035e4c4d49df4a06001b392c8014105`;
- FastDE `fastde abundance`, using the MIL bag encoder and gated attention
  pooling implemented in `src/fastde/abundance_mil.py`.

The synthetic reference fixture encodes cell type as one-hot single-cell
features. Original scSurvival receives these cells through
`feature_flavor='Custom'`; FastDE receives the same biological signal through
sample bags built from `sample x cell_type` counts.

This benchmark checks interface compatibility, top signal recovery, survival
C-index, risk-score agreement, wall time, and peak RSS. It is not a claim that
abundance-only bags reproduce scSurvival's expression AE/PCA workflows.

Command:

```bash
PYTHONPATH=src python scripts/fastde/benchmark_abundance_scsurvival_reference.py \
  --output-dir tmp/fastde_abundance_scsurvival_reference \
  --fastde-python /mnt/sdb/xzh/Vproject/SCOOP/.venv-scoop-fast/bin/python \
  --scsurvival-python /mnt/sdb/xzh/Vproject/SCOOP/.venv-scoop-omicverse/bin/python \
  --scsurvival-source /tmp/scSurvival_ref \
  --force-cpu \
  --epochs 60 \
  --learning-rate 0.02 \
  --n-samples 24 \
  --cells-per-sample 120 \
  --n-celltypes 4
```

Outputs:

- `benchmark_summary.tsv`
- `benchmark_consistency.json`
- `fastde/fastde_reference_metrics.json`
- `scsurvival/scsurvival_reference_metrics.json`
- framework-specific prediction tables

## Stress Test

The active stress entry point is:

```bash
PYTHONPATH=src python scripts/fastde/stress_abundance_mil.py \
  --output-dir tmp/fastde_abundance_mil_stress \
  --n-samples 72 \
  --n-celltypes 8 \
  --cells-per-sample 1000 \
  --max-instances-per-sample 256 \
  --epochs 40 \
  --learning-rate 0.02
```

It runs binary, multiclass, and survival modes end to end, records wall time and
peak RSS, and writes:

- `stress_summary.tsv`
- `stress_summary.json`
- one FastDE abundance output directory per mode

Local run on 2026-06-16:

```bash
PYTHONPATH=src .venv-scoop-fast/bin/python \
  scripts/fastde/stress_abundance_mil.py \
  --output-dir tmp/fastde_abundance_mil_stress \
  --n-samples 72 \
  --n-celltypes 8 \
  --cells-per-sample 1000 \
  --max-instances-per-sample 256 \
  --epochs 40 \
  --learning-rate 0.02
```

| Mode | Wall time | Peak RSS | Top cell type | Key metric |
| --- | ---: | ---: | --- | ---: |
| binary | 4.81s | 1.39 GB | ResponderCells | AUC 1.000 |
| multiclass | 3.33s | 1.40 GB | Other0 | macro AUC 1.000 |
| survival | 4.06s | 1.69 GB | RiskCells | C-index 0.781 |

Fixture details:

- 72 samples;
- 8 cell types;
- 1,000 cells per sample;
- bags capped at 256 sampled instances per sample;
- binary response driven by `ResponderCells`;
- survival risk driven by `RiskCells`;
- subtype B driven by `Other0`, so `Other0` is the expected top multiclass
  signal in this fixture.

The old linear benchmark numbers are intentionally not carried forward in this
document. They measured a different method and are no longer a valid description
of the default FastDE abundance implementation.
