# Example One: scSurvival for DAP12 / GSE274546

Search terms: scSurvival / OS / Cox / MIL / attention / cellular risk / myeloid / MYD88.
Verification date: 2026-10-07. Actual existing outputs were read here, without retraining.

## Questions to answer

Which expression states within myeloid cells carry patient survival risk? Then examine these states' identities,
MYD88 expression, and per-patient capture proportions to determine whether a “risk state” and “more cells in a population” point in the same direction.
This is expression-feature-driven scSurvival; proportion regression cannot stand in for it.

## Actual inputs and execution

Project root: `/mnt/sdb/xzh/Vproject/DAP12`.

| Object | Actual path within the project | Contents |
| --- | --- | --- |
| Cell features | `result/GSE274546_myeloid_survival/scsurvival_input/gse274546_myeloid_cell_features.csv.gz` | barcode, sample, patient, myeloid state, MYD88, Harmony features |
| Clinical table | `result/GSE274546_myeloid_survival/os_primary_paired_samples.csv` | 56 samples/patients; OStime, OS; 46 events |
| Original execution script | `script/run_gse274546_scsurvival_myd88.py` | Directly calls the original package's `scSurvivalRun` |
| Cell output | `result/GSE274546_myeloid_survival/scsurvival/cell_scsurvival_scores_primary_os.csv` | 39,842 cells, 54 patients actually covered |
| Patient output | `result/GSE274546_myeloid_survival/scsurvival/sample_scsurvival_hazards_primary_os.csv` | 54 rows, 45 events; time/status/patient_hazards |
| State output | `result/GSE274546_myeloid_survival/scsurvival/cluster_scsurvival_myd88_summary.csv` | Risk, attention, and expression summaries for 15 clusters |

Actual outputs lack `P28T1` and `P75T1` from the clinical table; do not report 56 as the number of patients covered by model outputs.
This check only verifies the coverage difference and does not infer why they are missing.

The old script selects the first 30 Harmony dimensions and places them in `obsm['X_harmony']` after StandardScaler;
it uses `feature_flavor='Custom'` and `feature_key='X_harmony'`. The single all-zero column in AnnData
X is an interface placeholder; actual inputs are expression-derived Harmony features, so do not claim that the model was trained on zero expression.
Parameters include 5 heads, dropout 0.3, epochs 240, patience 60, lr 0.001,
instance_batch_size 2000, sample_balance=True, predict_nMC=32, and validate=False.
These records explain old artifacts, rather than prescribing fixed optimal parameters for other cohorts.

## What the existing numbers look like

Read the [complete 15-cluster table](data/dap12-scsurvival-summary.csv) directly.
Three rows are shown below; risk is the pooled-cell mean from the old table, not a clinical HR.

| State | Cells | Samples detected | mean hazard_adj | mean attention |
| --- | ---: | ---: | ---: | ---: |
| TAM_SELENOP_MRC1 | 2,779 | 45 | −2.5290 | 0.6277 |
| Myeloid_Inflammatory_CCL3_IL1A (cluster 5) | 2,522 | 51 | −2.3242 | 0.7479 |
| TAM_APOC2_MT3 | 150 | 30 | −5.0002 | 0.7982 |

The old state table additionally joins results from a separate set of proportion models: the SELENOP/MRC1 row has OS HR=0.7435,
P=0.0430; the CCL3/IL1A row has OS HR=0.7782, P=0.2620, but its
OR_no_later_recur=0.4614, P=0.0290. These come from
`result/GSE274546_myeloid_survival/myeloid_scSurvival_like_results.csv`, rather than
scSurvival model outputs; joining the tables alone also does not establish that they use the same patients and effect units.
For the next comparison, verify coverage and units against each model's inputs/receipt.

## Transferable analytical and graphical moves

Begin with a continuous field of cellular risk, follow with “patient × state” risk distributions, attention, and identity/function
membership profiles, and finally compare state proportions and clinical effects within the same patient set. Risk distributions, expression identity, and
captured quantities each address a different question; establish a visual focus with shared state IDs and ordering.
Narrow clinical tracks and proportion effects can sit beside the state matrix; original tissue images provide spatial localization, with experimental slots reserved for
sorting/localization and functional perturbation. Do not draw attention as another undefined “cell proportion.”

This example yields three direct lessons:

- Low OS risk and less recurrence do not necessarily point in the same direction; read endpoints separately.
- APOC2_MT3 has very low risk, but its top markers are predominantly mitochondrial genes; interpret it together with quality and identity,
  rather than automatically naming it “protective macrophages” because it has the lowest risk.
- Old outputs with `validate=False` are fitted risk profiles. If a new project claims predictive ability, add
  patient-level split held-out risks and clinical baselines; do not treat old fitted values as external validation.

The old script also contains group-selection rules based on risk quantiles, proportion direction, and IFN exceptions; these are historical screening records,
and must not become Scoop's automatic biological annotation rules. Reuse the organizational sequence of cell risk→state function→patient comparison,
while establishing state interpretation and evaluation anew for the current project.

See the [case manifest](data/source-manifest.json) for sources and hashes. The algorithm reference is the
[original scSurvival repository](https://github.com/cliffren/scSurvival).
