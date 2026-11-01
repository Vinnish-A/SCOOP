---
name: scoop
description: Run and review single-cell and spatial omics with SCOOP, including QC, evidence-based annotation, CNV, programmes, atlas integration, spatial deconvolution and supporting multimodal data alignment. Also use for SCOOP runtime maintenance.
---

# SCOOP

Computation supplies evidence; the analyst/reasoning model adjudicates biology.
Read [measurement and evidence principles](references/docs/design-principles.md)
for analytical decisions, and only the relevant stage reference below.

Complete the provenance-linked metadata and dataset mapping in the
[metadata and batch contract](references/docs/design-principles.md#元信息与批次契约)
before integration. Latent correction, expression correction and statistical
batch control are separate recorded operations.

For a project-wide preprocessing/alignment request, first inventory all available
modalities and original measurement units using [multimodal preparation](references/stages/multimodal.md).
Look beyond the current RNA/spatial manifest for omitted source assets. Reuse
verified work and process executable gaps; a consistent ID bridge alone does not
complete native QC, biological annotation or downstream analysis.

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
`scoop status --run-root RUN --verify` distinguishes computation, validation,
interpretation and coverage. Register verified project-adapter artifacts through
`scoop register-result`; file presence or a successful exit is not biological approval.

## Routes

| Task | Read |
| --- | --- |
| QC / inherited lineage QC | [QC](references/stages/qc.md); [assay ranges](references/docs/assay-preprocessing.md) for imaging/spot assays |
| Embedding / parallel clustering | [core](references/stages/core.md) |
| Programme discovery | [programmes](references/stages/programs.md) |
| Cell identity / subtype review | [annotation workflow](references/markers/annotation-workflow.md), [naming](references/markers/naming-strategy.md), then tissue/lineage references |
| Malignant/normal evidence | [tumor rules](references/markers/tumor-general.md) |
| Multi-study scVI / approved scPoli | [atlas](references/docs/large-atlas-integration.md) |
| Multimodal inventory; ATAC/ChIP, ADT/IMC, bulk/protein, perturbation and clinical bridges | [multimodal preparation](references/stages/multimodal.md) |
| Spatial deconvolution | [spatial](references/stages/spatial.md) |
| Communication candidates | [CCC](references/stages/ccc.md) |
| Condition DE | [DE](references/stages/differential-expression.md) |
| Proportions / scores / survival / response | [patient outcomes](references/stages/outcome-association.md) |

For spatial deconvolution, cover **all reviewed cell subgroups across the whole
tissue** in the fitted reference, not only the focal lineage. Fit the fine
reference first; aggregate its weights for broad displays when needed. Keep
unresolved populations explicit and report unassessed subgroup gaps; do not
replace reviewed subtypes with broad classes for convenience. See the spatial
route for the reference coverage contract.

## VEngine中的SCOOP

SCOOP是V体系的基础层：清点、测量、预处理、跨数据集与跨模态对齐、已审核身份与回执在此完成；VEngine的发表分析（CH方法链）在这些对象上计算，其索引列出每条链所需的SCOOP路线；VStory与OmniFigure承担论证与成图。在这里完成或核实方法链所需的路线，再交付对象、键与状态；方法链的结果表用`scoop register-result`登记；SCOOP尚无的阶段记录为缺口，被多个项目复用后提升为SCOOP阶段。SCOOP不依赖这些技能，上述内容在没有它们时同样适用。本指南随每个运行时发布（`${SCOOP_HOME:-~/.scoop}/current/skill/zh/`），经VEngine读取，不再作为独立的agent技能安装。

## Biological decisions

The main agent coordinates delegated execution and evidence collection and
personally reads the actual clustering, DEG, UMAP distances/neighbourhoods and
tissue markers to adjudicate biologically meaningful broad cell classes. Numeric
clusters are not themselves biological classes; reviewed classes may merge or
split them. Credible clustering, conventional DEG and tissue markers can establish
a broad class without NMF. The main agent personally reviews the actual parent
NMF and all other required evidence before formal fine-subtype adjudication; a
subagent's proposed label is not a final decision. See the
[annotation workflow](references/markers/annotation-workflow.md) for the full
readiness and evidence gates. SingleR,
CellTypist, scANVI or label transfer must not assign final identities. Scores and
programmes provide evidence. Keep identity, state, malignancy and confidence
separate; retain unresolved labels. Validate structured decisions and commit
through the deterministic annotation interface.

Missing evidence is `evidence_pending` workflow progress, not an unresolved
biological identity: actively retrieve source data, tissue markers and literature
or run necessary analyses, documenting unavailable measurements or execution
limits and continuing feasible work. Do not submit a final type while required
evidence is missing. Reserve identity `unresolved` for adequate evidence already
reviewed that remains conflicting or mixed. This workflow status is not a new
runtime annotation-schema enum.

CNV-low does not establish normality; normal references require reviewed eligible
identities. scVI learns representations. scPoli requires explicit current approval
bound by SHA-256 to the report, reviewed labels, input and labeled indices; only
high-confidence reviewed cells supervise it. Patient outcomes never decide identity.
