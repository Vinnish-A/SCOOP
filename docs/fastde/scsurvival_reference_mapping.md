# scSurvival Reference Mapping

## Source Located

- Repository: https://github.com/cliffren/scSurvival
- Local inspection path: `/tmp/scSurvival_ref`
- Commit inspected: `a76de0a00035e4c4d49df4a06001b392c8014105`
- Package version: `1.3.0` from `setup.py`
- License: GPL-3.0 from `setup.py` and `LICENSE`

The source was not present in the SCOOP repository, submodules, or the existing
Fast and OmicVerse virtual environments. It was found on public GitHub and
cloned for inspection.

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

## Copied Or Ported

No GPL source file was copied verbatim into SCOOP. The implementation ports the
public API shape and survival-loss behavior at a clean-room level:

- Cox partial likelihood with right censoring and stable log-risk-set sums.
- scSurvival-compatible names in `src/fastde/scsurvival_compat.py`:
  `ScSurvivalDataset`, `ScSurvivalModel`, `ScSurvivalTrainer`,
  `ScSurvivalResult`.
- Similar runner concept through `fastde abundance`.

## Reimplemented

The FastDE abundance module reimplements the method around SCOOP's requested
input: a sample-by-cell-type abundance matrix. This is not the same data model
as original scSurvival's patient-by-single-cell expression bags.

Reimplemented files:

- `src/fastde/abundance_data.py`: H5AD/sample count aggregation.
- `src/fastde/abundance_design.py`: abundance transforms and design matrices.
- `src/fastde/abundance_loss.py`: Cox, classification, and multinomial losses.
- `src/fastde/abundance_model.py`: deterministic NumPy/sklearn model backend.
- `src/fastde/abundance_train.py`: training wrapper.
- `src/fastde/abundance.py`: mode-specific results, predictions, metrics, and
  manifests.
- `src/fastde/abundance_cli.py`: `fastde abundance` CLI.

## Differences From Original scSurvival

- Input is sample-level cell-type or cell-state counts/proportions, not
  single-cell expression bags.
- Default transform is centered log-ratio with pseudocount.
- The default backend is NumPy/SciPy/sklearn so it runs in SCOOP's Fast
  environment without requiring PyTorch or GPU.
- Survival mode uses a linear Cox risk model on transformed abundance features
  and optional covariates.
- Binary/multiclass/continuous modes are added for FastDE abundance and are not
  part of the inspected scSurvival public survival runner.
- Large per-sample/cell-type tables are written as external TSV/JSON artifacts
  and are not stored inside H5AD.
- The auxiliary abundance reconstruction loss is documented in the manifest but
  not enabled in the current deterministic NumPy backend.

## Compatibility Gap

The implementation is scSurvival-like, not a direct scSurvival clone. The main
reason is methodological: FastDE abundance tests sample-level proportions,
whereas original scSurvival learns from raw single-cell expression for each
patient. The shared contract is survival/outcome association from cohort data,
Cox-style survival loss, comparable diagnostics, and compatibility class names.
