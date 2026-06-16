# FastDE Abundance vs scSurvival Reference Benchmark

This benchmark compares `fastde abundance --mode survival` with the original
public scSurvival implementation on a controlled synthetic cohort.

The comparison fixture encodes cell type as one-hot single-cell features:

- scSurvival receives per-cell custom features through `feature_flavor='Custom'`.
- FastDE abundance receives the same cells aggregated to sample-by-celltype
  counts/proportions.

This is the closest fair comparison for the new FastDE abundance scope. It does
not claim equivalence to full original scSurvival AE/PCA expression modeling,
because FastDE abundance intentionally tests sample-level proportions instead of
single-cell expression bags.

Command:

```bash
PYTHONPATH=src python scripts/fastde/benchmark_abundance_scsurvival_reference.py \
  --output-dir tmp/fastde_abundance_scsurvival_reference \
  --scsurvival-python /mnt/sdb/xzh/Vproject/SCOOP/.venv-scoop-omicverse/bin/python \
  --scsurvival-source /tmp/scSurvival_ref \
  --force-cpu
```

Reference environment note: original scSurvival requires `torch`, `scanpy`,
`scikit-learn`, and `lifelines`. In the local run, `lifelines==0.30.3` was
installed into the separate OmicVerse/reference environment, not the Fast
environment.

Outputs:

- `benchmark_summary.tsv`
- `benchmark_consistency.json`
- `fastde/fastde_reference_metrics.json`
- `scsurvival/scsurvival_reference_metrics.json`
- framework-specific prediction tables

Metrics:

- wall time;
- peak RSS;
- survival C-index;
- top cell-type signal;
- Spearman correlation between FastDE and scSurvival sample risk scores;
- speedup and RSS ratio.

## Local Synthetic Result

Run date: 2026-06-16.

Command:

```bash
PYTHONPATH=src .venv-scoop-fast/bin/python \
  scripts/fastde/benchmark_abundance_scsurvival_reference.py \
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

Fixture:

- 24 samples;
- 2,880 cells;
- 4 cell types;
- survival risk driven by `RiskCells` abundance;
- CPU-only for both frameworks;
- scSurvival used `feature_flavor='Custom'` with one-hot cell-type features.

| Framework | Wall time | Peak RSS | C-index | Top cell type |
| --- | ---: | ---: | ---: | --- |
| FastDE abundance | 0.036s | 186.2 MB | 0.882 | RiskCells |
| scSurvival reference | 15.252s | 942.8 MB | 0.878 | RiskCells |

Consistency:

- sample risk Spearman: `0.869` (`p=3.71e-08`);
- top cell-type signal match: `true`;
- speedup vs scSurvival: `424.8x`;
- scSurvival/FastDE peak RSS ratio: `5.06x`.

Interpretation:

FastDE abundance is much faster here because it solves the intended sample-level
proportion problem after aggregation. Original scSurvival carries single-cell
MIL/attention training overhead even when the supplied custom features are just
cell-type indicators. The two methods are consistent on the designed signal in
this benchmark, but this should not be interpreted as full equivalence to
scSurvival's expression-level AE/PCA workflows.
