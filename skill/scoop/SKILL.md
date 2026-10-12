---
name: scoop
description: Execute and review auditable single-cell and spatial omics workflows with SCOOP. Use for H5AD validation, per-sample QC and doublet detection, embedding and clustering, FastCNMF programmes, evidence-based annotation, tumor CNV review, large-atlas scVI annotation and scPoli gating, spatial deconvolution, cell-cell communication, pseudobulk differential expression, abundance and patient-outcome analysis (scSurvival, relative proportions, UCell, Cox or response), delivery, or maintenance of the SCOOP runtime and references.
---

# SCOOP

Keep biological judgment in the Skill and deterministic computation in the
runtime. Read [first principles](references/docs/design-principles.md) before
analysis: preserve measurement meaning; bypass unusable components and record
what actually ran. Defaults and technical pass flags do not establish biology.

## Execution boundary

Prefer `${SCOOP_HOME:-$HOME/.scoop}/bin/scoop` and its managed environment when
available. SCOOP Home must resolve to an absolute path. `doctor --json` discovers
dependencies; test the actual backend on a representative input before expensive
work. Inspect coverage, batch/epochs/warmup, memory and scheduling.

If the interface is missing, hangs or scales poorly, use the real API or a
project-local adapter, with isolated dependencies when needed. Follow the
[runtime contract](references/runtime.md); record why, actual method/environment,
inputs/parameters, outputs and checks. Do not silently switch algorithms or
modify the shared release to repair one project. A bypass is not managed-stage
completion or biological approval.

- Keep raw data read-only, integer counts in `layers["counts"]`, project outputs
  outside SCOOP Home, large tables external and stable state/registry in H5AD.
- Reuse validated complete artifacts and checkpoints; file existence is not
  completion. Invalidate affected descendants only. Save/interface failures do
  not require rerunning valid analysis.
- Validate evidence-backed structured annotation, then commit deterministically;
  a project committer may replace an unusable CLI, not the biological review.

Use [workflow contract](references/workflow-contract.md) for dependencies,
receipts and scoped review; [project structure](references/docs/project-structure.md)
when creating a run. Use `scoop plan --workflow <workflow> --config <run.yaml>`
when supported, otherwise record the project execution plan.

## Route the workflow

Choose the narrowest useful route. Workflow names do not guarantee installed
command support; missing adapters are handled under the execution boundary.

| Request | Route | Read before execution |
| --- | --- | --- |
| QC, embedding, programmes, annotation | `standard` | [modules](references/docs/modules-and-outputs.md), relevant [QC](references/stages/qc.md), [core](references/stages/core.md), [programmes](references/stages/programs.md), [annotation](references/stages/annotation.md) |
| Malignant/normal review | `tumor` | [tumor rules](references/markers/tumor-general.md) and matching tissue reference |
| Multi-study atlas, scVI/scPoli | `atlas` | [large atlas](references/docs/large-atlas-integration.md), [annotation workflow](references/markers/annotation-workflow.md) |
| Spatial deconvolution/graph | `spatial` | [spatial](references/stages/spatial.md), [parameters](references/docs/parameter-tuning.md) |
| Candidate communication | `ccc` | [CCC](references/stages/ccc.md) |
| Condition DE | `differential-expression` | [DE](references/stages/differential-expression.md) |
| Proportions, UCell, scSurvival, survival/response | `differential-expression` | [outcomes and real examples](references/stages/outcome-association.md) |

## Annotation and supervision

Read [annotation workflow](references/markers/annotation-workflow.md) and
[naming](references/markers/naming-strategy.md), then relevant
[major](references/markers/major-cell-types.md), [T](references/markers/t-cells.md),
[myeloid](references/markers/myeloid.md), [fibroblast](references/markers/fibroblast.md),
[BRCA](references/markers/tumor-brca.md), [CRC](references/markers/tumor-crc.md) or
[GBM](references/markers/tumor-gbm.md) references.

The analyst/reasoning model adjudicates identity using tissue-matched atlases,
verifiable literature and actual marker, sample and cluster evidence. SingleR,
CellTypist, scANVI and automatic label transfer must not determine final labels.
Keep identity, state, malignancy and confidence separate; outcomes never decide
identity. CNV-low is not normal by default; reference identities require review.

scVI learns representations, not labels. Its report must show subtype/neighbour
structure, markers, reference conflicts, programmes, QC/confidence and proposed
scPoli supervision. scPoli requires a current approved gate binding the report,
labels, input H5AD and labeled indices by SHA-256. Exclude unreviewed legacy/scANVI
labels and medium/low-confidence, ambiguous or suspicious cells from supervision.
Tool bypass cannot substitute for this approval.

## Delivery

Return the biological result, scope and material gaps through one clean current
entry with validated objects/source tables, actual runtime/adapter, review state
and reproduction paths. Preserve cell/gene alignment, counts, clinical joins and
embedding/graph descriptors; validate storage without repeating the analysis.
