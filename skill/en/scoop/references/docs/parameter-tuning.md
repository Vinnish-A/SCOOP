# 05. Parameter Tuning Rules

The principle of parameter tuning is to first determine whether a problem comes from the data, the model, or biology, then change as few parameters as necessary.

## QC

Do not use globally fixed thresholds. Use each sample's own distribution. Only the following situations trigger manual review:

- More than 30% of a sample fails QC;
- Ribosomal-high exceeds 20% and affects PCA;
- High mt and low genes do not co-occur, yet a single metric produces a fail call;
- A condition almost disappears because of QC;
- A spatial FOV/tile coincides completely with QC failure.

## Batch correction

Harmony 2.0 is the default. Correct only technical keys such as `sample_id`, `library_id`, `batch_id`, and `chemistry`. Do not correct `condition`, `disease`, `treatment`, or `lineage`.

If batch and condition are highly confounded, do not perform correction; handle this in downstream DE or stratified analysis instead.

## Leiden resolution

Do not run only one resolution. The default search range is `0.25–1.5`, with 5 seeds. The default strategy is coarse-to-fine:

1. Start with a coarse sweep over `[0.25, 0.75, 1.25, 1.5]` and 2 seeds;
2. Select the lowest usable resolution that is stable across seeds and produces neither too few nor excessively fragmented clusters;
3. Run a fine sweep within `±0.25` of that resolution in `0.125` steps, using all 5 seeds;
4. Ultimately select the lowest stable resolution, rather than the resolution with the most clusters.

When selecting, examine:

- Marker support;
- Graph coherence;
- NMF programme coherence;
- Seed stability;
- Absence of dominance by stress/ribo/cell-cycle effects.

## Program NMF

Use the FastCNMF profile by default: `programs.method=fastcnmf`, `max_iter=50`, 20 seeds, and a K sweep of `[5, 8, 10, 12, 15, 20]`. `max_iter=50` is the more reliable default in the current internal S2 and spatial S1 benchmarks; if programmes are unstable, first stratify by broad lineage and rerun, rather than immediately increasing K or the iteration count substantially.

## Stress / ribosomal / proliferation

Do not remove these cells by default. First calculate within-sample z-score/MAD outliers, then examine whether they dominate PCs, cluster across lineages, have clear parent identities, or co-occur with doublets/ambient/low counts/high mt/bad FOVs.

Decision order:

1. Retain the biology embedding;
2. If cell-type clustering is affected, exclude the relevant genes from identity HVGs;
3. Or regress only the identity embedding;
4. If the parent identity is clear, retain the signal as cell_state;
5. Remove cells only when multiple lines of evidence establish an artifact.

## RCTD-py mode

The default for low-resolution spots/ROIs/large bins is `full`, because it estimates a continuous mixture of all cell types. `multi` is for sparse interpretation and sensitivity analysis, particularly when each spatial unit is expected to contain only a few dominant cell types.

## CCC

FastCCC is a screen, not proof of mechanism. Complex receptors, integrins, ECM, MHC, cytokine receptors, and similar interactions must undergo complex-aware validation. Spatial projects must also check whether sender and receiver are adjacent or co-localized.
