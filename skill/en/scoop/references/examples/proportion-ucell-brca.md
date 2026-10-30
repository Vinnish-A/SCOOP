# Example Two: Relative Proportions, UCell, and pCR in BRCA / ARTEMIS

Search terms: zzu_gp / relative fraction / parent lineage / UCell / rare state /
pCR / OR / LAM / LPL / forest / clinical annotation.
Verification date: 2026-10-07. Existing patient tables and regression results were reused, with no new fitting.

## Questions to answer

How do a state's “fraction of its parent population” and the “strength” of its expression programme in a patient each relate to
neoadjuvant treatment response? Measuring quantity and function separately allows LPL to emerge naturally from the ecological context of LAMs.

## Actual cohort and two measurements

Project root: `/mnt/nfs_data/xzh/BRCA`.
After QC: 417,698 cells and 101 patients; the original 49-state forest plot included 82 patients with available outcomes,
45 pCR and 37 RD. Freeze original states and reviewed identities; do not obtain a clinical direction by changing labels.

The earlier method comes from `zzu_gp/R/0Fig1.5.R:692–780` and is implemented in this project as follows:

| Measurement | Actual definition | Must not be conflated with |
| --- | --- | --- |
| Common-state proportion | Per-patient state_count / corresponding parent-lineage count | Fraction of all cells in the whole sample |
| Rare-state UCell | Fixed published markers; all-gene ranking; mean of scores strictly >Q70 among all cells of each patient | Positive-cell fraction, counts pseudobulk, actual LAM abundance |
| LAM count comparator | Original C1Q/CD36 state ∩ approved Myeloid, divided by all original mac-* cells in the same scope | UCell programme strength across all patient cells |
| Myeloid LPL output | log2(1+CPM) of counts pseudobulk within approved Myeloid | LAM-specific expression, UCell, cell proportion |

The earlier rare-state rule was “mean whole-cell fraction among patients in whom that state is detected <0.015,” not a zero-filled mean across the entire cohort.
C1Q=0.013348 and CD36=0.007004, so both rows in the earlier forest plot use UCell.
This is a record of the project method; not every new project must adopt this cutoff.

Actual UCell ranks all 32,354 genes, with maxRank=1500, average ties, and missing-gene
imputation; the 49 fixed signatures contain 14–20 published markers, with 20 in each LAM signature. The source is the counts layer of canonical
`yan_counts.h5ad`, rather than an integrated matrix containing only HVGs. The main summary uses all cells from each patient,
with within-parent-lineage upper tails saved separately. Validation against R UCell 2.14.0 has been performed on 72 actual cells ×49 signatures, with a maximum
absolute error of 2.220446e−16. This checks equivalence of a historical implementation; it does not imply equivalence across different markers/scopes.

Common-state coverage depends on whether the parent population is present: fill in 0 when the parent is present but the target is absent; the proportion is NA when the parent is absent.
Thus, patient counts for common states in the earlier forest plot include T=82, Mye=79, B=70, and Fibro=20; they cannot all be reported as 82.

## Existing clinical results

The forest plot consistently uses pCR=1/RD=0, with continuous exposures scaled by 1 SD within the model's patient set.
If an original output treats RD as the event, invert it to obtain OR_pCR and exchange the CI bounds, leaving P unchanged.

| Actual metric | Patients | Model | OR_pCR | 95% CI | nominal P |
| --- | ---: | --- | ---: | --- | ---: |
| C1Q UCell upper-tail mean across all patient cells | 82 | unadjusted | 0.7662 | 0.4921–1.1931 | 0.2385 |
| CD36 UCell upper-tail mean across all patient cells | 82 | unadjusted | 0.8587 | 0.5546–1.3296 | 0.4947 |
| C1Q / Mac raw proportion | 70 | unadjusted | 1.7458 | 1.0415–2.9265 | 0.03449 |
| C1Q / Mac raw proportion | 70 | composition | 1.6866 | 0.9774–2.9101 | 0.06038 |
| Myeloid LPL log2(1+CPM) | 70 | composition | 0.5025 | 0.2902–0.8702 | 0.01404 |

UCell in 82 patients and proportions in 70 patients do not have the same exposure/coverage. Comparisons for the latter two measurement types within the same 70 patients include
37 pCR and 33 RD; the full set defining the parent population comprises 22,744 eligible macrophages, with a parent-population requirement of ≥30 per patient.
The rare-LAM numerator itself has no ≥30 cutoff; CD36 is nonzero in only 8 patients, with 62 true zeros.
Composition covariates are fixed as the Myeloid, T, and Epithelial fractions of all cells.

The conclusion actually supported is that C1Q capture proportion and myeloid LPL output have different clinical directions; in these data,
the response association is clearer for LPL output. UCell did not show a significant response association for overall LAM programme strength,
and must not be relabeled as a significant LAM proportion. Fixed-logit sensitivity analysis is another measurement; its P
cannot replace the P for the raw proportion. Another model on the same patients is also not independent cohort validation.

## Directly usable tables and sources

- [Main forest table for all 49 states](data/brca-state-primary-or.tsv): includes original state IDs, measurement type,
  OR/CI/P, coverage, original RD direction, and q; it retains more than just the two LAM rows.
- [Quantity/function comparison within the same 70 patients](data/brca-lam-lpl-comparison.tsv): all 10 original comparison rows.
- [LAM rare-state method audit](data/brca-lam-rare-audit.tsv): actual cell counts, proportions, and signature sources.
- Original 49-state results: `results/fig1_vstory_20261005/rebuild_scrna/`, where
  `forest49_primary_OR_pCR.tsv`, `analysis_config.json`, and
  `REBUILD_SCRNA_METHODS.zh.md` preserve main results and parameters.
- Original patient count table: `results/LAM_outcome_EDA_20261005/SC/patient_fixed_LAM_allmac_estimands.tsv`;
  original model tables and reports share that directory, and code is in `analysis/LAM_outcome_EDA_20261005/SC/`.

See the [case manifest](data/source-manifest.json) for numeric copies and source paths/hashes.

## Building an information-dense main figure

Use a state identity/function heatmap and a clinical forest plot with a shared index, semantically group common proportions and rare scores,
and add narrow patient tracks and clinical-metric comparisons. Introduce LPL through lipid-programme membership and lineage expression before the forest plot,
and use a same-patient quantity–function comparison at the side to explain the focus on LPL; do not switch state names merely for significance.
Connect single-cell outcomes to native tissue, niches, and within-patient spatial comparisons using actual paired crosswalks.
Align the actual data edges of the heatmap and forest plot, share patient order across clinical tracks, and fit legends into gaps.
Reserve experimental slots for LAM/LPL localization and functional validation; do not present computed scores as experimental readouts.

For a new project, rebuild the patient×state table from frozen identities and the complete gene table, then output OR or HR according to the actual endpoint.
The earlier 0.015/Q70/top20 parameters are reviewable examples; after changing the parent population, markers, or
summary, give the new measurement its own name and complete results, rather than claiming that only the computational backend of the same model was changed.
Algorithm reference: [official UCell repository](https://github.com/carmonalab/UCell).
