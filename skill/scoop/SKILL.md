---
name: scoop
description: Run and review single-cell and spatial omics with SCOOP, including QC, clustering, programmes, evidence-based annotation, tumor CNV, atlas integration, spatial deconvolution and patient-level associations. Also use for SCOOP runtime maintenance.
---

# SCOOP

Computation supplies evidence; the analyst/reasoning model adjudicates biology.
Read [measurement and evidence principles](references/docs/design-principles.md)
for analytical decisions, and only the relevant stage reference below.

## Execution

Use `${SCOOP_HOME:-$HOME/.scoop}/bin/scoop` or `bin/python`. The CLI and API share
one locked dependency environment; application/Skill updates reuse it. Read
[runtime](references/runtime.md) for installation, capabilities and rollback.

`doctor --json` discovers dependencies/adapters, not successful analysis. Test the
actual backend on representative data before costly work. Check assay, library
scope, batch/epochs/KL warmup, memory and task-level concurrency. For large H5AD,
use [selective reads and inherited QC](references/docs/data-boundaries.md).

If an interface is missing, hangs or scales poorly, use the real API or a project
adapter and record method/environment, reason, inputs/parameters, outputs and
checks. Isolate genuinely conflicting dependencies; do not modify a shared release
or silently switch algorithms. A bypass does not grant biological approval.

Preserve raw counts, cell/gene order, clinical joins and source files. Resume
validated checkpoints; a save failure does not require refitting valid models.
Keep large tables external and preserve executable graph/embedding metadata.
[Workflow receipts](references/workflow-contract.md) describe actual execution
and review state; use the narrowest route that answers the question.

## Routes

| Task | Read |
| --- | --- |
| QC / inherited lineage QC | [QC](references/stages/qc.md); [assay ranges](references/docs/assay-preprocessing.md) for imaging/spot assays |
| Embedding / parallel clustering | [core](references/stages/core.md) |
| Programme discovery | [programmes](references/stages/programs.md) |
| Cell identity / subtype review | [annotation workflow](references/markers/annotation-workflow.md), [naming](references/markers/naming-strategy.md), then tissue/lineage references |
| Malignant/normal evidence | [tumor rules](references/markers/tumor-general.md) |
| Multi-study scVI / approved scPoli | [atlas](references/docs/large-atlas-integration.md) |
| Spatial deconvolution | [spatial](references/stages/spatial.md) |
| Communication candidates | [CCC](references/stages/ccc.md) |
| Condition DE | [DE](references/stages/differential-expression.md) |
| Proportions / scores / survival / response | [patient outcomes](references/stages/outcome-association.md) |

## Biological decisions

Review DEG, positive/exclusion markers, cluster structure, sample/QC distribution
and relevant literature; propose broad identities before subtypes. SingleR,
CellTypist, scANVI or label transfer must not assign final identities. Scores and
programmes provide evidence. Keep identity, state, malignancy and confidence
separate; retain unresolved labels. Validate structured decisions and commit
through the deterministic annotation interface.

CNV-low does not establish normality; normal references require reviewed eligible
identities. scVI learns representations. scPoli requires explicit current approval
bound by SHA-256 to the report, reviewed labels, input and labeled indices; only
high-confidence reviewed cells supervise it. Patient outcomes never decide identity.
