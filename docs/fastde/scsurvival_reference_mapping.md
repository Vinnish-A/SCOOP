# scSurvival Reference Mapping

## Source Located

- Repository: https://github.com/cliffren/scSurvival
- Local inspection path: `/tmp/scSurvival_ref`
- Commit inspected: `a76de0a00035e4c4d49df4a06001b392c8014105`
- Package version: `1.3.0` from `setup.py`
- License: GPL-3.0 from `setup.py` and `LICENSE`

The source was not present in the SCOOP repository at the start of this work.
It was found on public GitHub, inspected locally, and then pinned as a
third-party submodule at `third_party/scSurvival`.

## Files Inspected

- `setup.py`
- `LICENSE`
- `scSurvival/__init__.py`
- `scSurvival/base_module.py`
- `scSurvival/loss_func.py`
- `scSurvival/scsurvival.py`
- `scSurvival/scsurvival_core.py`
- `scSurvival/scsurvival_module.py`
- `scSurvival/utils.py`

## scSurvival Behavior Relevant To FastDE

Original scSurvival models survival outcomes from single-cell cohort data. It
groups cells by patient, learns cell-level representations with an AE/VAE-style
module, pools patient-level risk with attention-like components, and optimizes a
Cox partial likelihood. It also reports cell/sample risk profiles.

Important source functions/classes:

- `scSurvivalRun`: public high-level runner.
- `scSurvival`: PyTorch training wrapper for cohort survival modeling.
- `scSurvivalCellModelAE`, `scSurvivalCellModelVAE`, `scSurvivalCellModel`:
  single-cell representation modules.
- `HazrdModel`: hazard/risk head.
- `cox_loss_func`: Cox partial negative log likelihood on samples sorted by
  descending survival time.
- `c_index` and `conditional_cindex`: survival evaluation.

## Copied, Vendored, Or Ported

The original scSurvival repository is tracked as a third-party submodule for
source provenance and direct reference:

- `third_party/scSurvival`
- pinned commit: `a76de0a00035e4c4d49df4a06001b392c8014105`

No GPL source file is copied verbatim into `src/fastde`. The FastDE default
backend ports the architecture at a clean-room level:

- Cox partial likelihood with right censoring and stable log-risk-set sums.
- sample-as-bag multiple-instance learning.
- shared instance encoder.
- gated attention pooling.
- sample-level task head.
- scSurvival-compatible names in `src/fastde/scsurvival_compat.py`:
  `ScSurvivalDataset`, `ScSurvivalModel`, `ScSurvivalTrainer`,
  `ScSurvivalResult`.
- Similar runner concept through `fastde abundance`.

## Reimplemented

The FastDE abundance module reimplements the method around SCOOP's requested
input while preserving the original sample-as-bag idea. H5AD inputs use cells as
bag instances with cell-type or cell-state one-hot features. Direct count-matrix
inputs are expanded into equivalent bag instances for compatibility.

Reimplemented files:

- `src/fastde/abundance_data.py`: H5AD/sample count aggregation.
- `src/fastde/abundance_design.py`: abundance transforms and design matrices.
- `src/fastde/abundance_loss.py`: Cox, classification, and multinomial losses.
- `src/fastde/abundance_mil.py`: scSurvival-style MIL encoder, gated attention
  pooling, task heads, and survival loss variants.
- `src/fastde/abundance_model.py`: deterministic NumPy/sklearn model backend.
- `src/fastde/abundance_train.py`: training wrapper.
- `src/fastde/abundance.py`: mode-specific results, predictions, metrics, and
  manifests.
- `src/fastde/abundance_cli.py`: `fastde abundance` CLI.

## Differences From Original scSurvival

- Input focuses on cell-type or cell-state abundance rather than raw gene
  expression.
- The default backend is PyTorch `scsurvival_mil`; the older NumPy/SciPy linear
  model remains available as `--abundance-backend linear`.
- Survival mode supports the original Cox partial likelihood plus optional
  rank-loss variants.
- Binary/multiclass/continuous modes reuse the same scSurvival-style
  encoder-attention-bag architecture with task-specific heads. These modes are
  FastDE extensions and are not part of the inspected scSurvival public survival
  runner.
- Large per-sample/cell-type tables are written as external TSV/JSON artifacts
  and are not stored inside H5AD.

## Compatibility Gap

The implementation is scSurvival-style rather than a direct clone. The main
reason is task scope: FastDE abundance tests cell-type/cell-state abundance,
whereas original scSurvival learns from raw single-cell expression for each
patient. The shared contract is sample-as-bag modeling, gated attention pooling,
survival/outcome association from cohort data, Cox-style survival loss,
comparable diagnostics, and compatibility class names.
