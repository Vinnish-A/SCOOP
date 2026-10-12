# Single-cell states and patient outcomes

Read this reference for scSurvival, relative cell proportions, UCell, Cox, pCR,
response, or a clinical forest plot. These are outcome questions, not condition DE.
Two worked examples use actual project tables:

- [DAP12 / GSE274546: scSurvival](../examples/scsurvival-dap12.md).
- [BRCA / ARTEMIS: relative proportions and UCell](../examples/proportion-ucell-brca.md).

## Choose the measurement before the model

| Biological question | Patient measurement | Outcome model | Result meaning |
| --- | --- | --- | --- |
| Which cell expression states carry survival risk? | Cell expression-derived features within patient bags | Original scSurvival attention-based multiple-instance Cox model | Patient risk and cellular risk/attention profiles |
| Is a captured cell state more prevalent within its lineage? | State cell count / parent-lineage cell count, per patient | Cox for time/event; logistic regression for binary response | HR or OR for a specified increment in proportion |
| Is a rare state's programme stronger despite sparse capture? | Fixed-marker UCell scores, summarized per patient in a stated cell scope | Cox or logistic regression | Association with programme strength, not cell abundance |
| Does a lineage produce more of a particular gene? | Lineage pseudobulk expression with an explicit library denominator | Cox or logistic regression | Gene output association; distinct from state count and UCell |

The first two examples cover complementary routes. scSurvival is supervised by
survival; proportion/UCell summaries can be defined without outcomes and then
tested. Binary response has no censoring time and is not an input substitute for
the original survival model. A composition-only binary MIL adaptation is a
different method, even if described elsewhere as “scSurvival-style.”

## Patient and measurement contract

1. Join cell → sample/timepoint → patient → endpoint. Report unmatched IDs. Keep
   each patient on one side of any train/test split; repeated samples do not add
   independent survival observations. Choose baseline or another biologically
   defined timepoint before patient aggregation.
2. Freeze the reviewed cell identity, state markers and parent lineage. Keep
   identity, functional state and learned risk in separate fields. A risk ranking
   does not reannotate cell types.
3. For proportions, retain true zeros when the parent is adequately measured;
   parent absent means an undefined ratio. A coverage threshold belongs to the
   denominator; demanding many target cells removes rare-state zeros. State the
   model scale, e.g. one SD of the raw ratio versus a fixed logit transform.
4. For UCell, record gene universe, marker coverage, maxRank, ties/missing-gene
   handling, scoring scope and patient summary. All-cell tail mean, within-lineage
   mean and fraction above a fixed threshold answer different questions. Do not
   rank only the marker genes or replace missing genes silently. A constant score
   can have an empty strict upper tail; preserve an explicit missing measurement.
5. Choose Cox for documented time/event, logistic regression for pCR/response,
   and an ordinal model for ordered RCB. State event direction and units. A
   clinical HR/OR is not a scSurvival cell score or an attention value. Keep main
   and sensitivity measurements named separately.

The zzu_gp/BRCA observed-patient mean fraction <0.015 and within-patient strict
>Q70 summary are documented project choices, not universal defaults. Their exact
application and the alternative count measurement are shown in the BRCA example.

## scSurvival deliverables

Keep the patient clinical table, cell-to-bag mapping, expression/feature source,
model configuration and patient split together. Export:

- Per-cell `cell_id`, `patient_id`, state, `hazard`, `hazard_adj`, `attention`.
- Per-patient time/event, fitted risk, and held-out risk when actually produced.
- State-level distributions and patient-balanced summaries, with marker/function
  profiles. Label a pooled cell mean as such; large cell yields otherwise dominate.
- Patient-level evaluation and comparison to a clinical/proportion baseline on
  the same split, when evaluating prediction. A fit with `validate=False` supplies
  fitted profiles; it does not supply a held-out performance estimate.

## Proportion / score deliverables

Export patient measurements and a long model table. Retain:

`patient_id, state, parent_lineage, measurement_kind, n_target, n_parent,
score_scope, score_summary, endpoint, event_direction, effect_unit,
n_patients, effect, CI_low, CI_high, P, model, covariates, status`.

Separate observed abundance, programme scores and gene output in the table. Keep
the full tested state set and coverage in the report; select the main figure for
its biological argument without changing the numerical meaning of a row.

## Figure interface

Build a shared state order across identity/function heatmap, parent-relative
abundance or score tracks, and a dense clinical forest. Include relevant clinical
comparators and narrow outcome tracks. State-specific expression or programme
evidence introduces a focal gene before its clinical association. Patient/state
matrices and original tissue context add different evidence; repeating the same
association with several plot forms does not.

Use aligned data boundaries, readable labels and legends in module gaps. Reserve
experiment slots for identity/localization and perturbation validation. Source,
aggregation, coverage and model details belong in the report, not oversized legend
blocks. Keep effect axes distinguishable when showing Cox and response together.

## Runtime scope

This reference does not register a CLI command. Check actual support in the
`differential-expression` route; if the adapter is absent or unusable, implement
the resolved patient measurement/model through the real package API or a
project-local adapter under the [runtime bypass contract](../runtime.md).
Legacy scripts are provenance, not automatically validated implementations.
Verify the actual inputs and outputs; do not stop at a plan merely because the
managed CLI lacks the method.

Method references: [original scSurvival](https://github.com/cliffren/scSurvival),
[UCell](https://github.com/carmonalab/UCell).
