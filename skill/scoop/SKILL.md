---
name: scoop
description: Execute and review auditable single-cell and spatial omics workflows with SCOOP. Use for H5AD validation, per-sample QC and doublet detection, embedding and clustering, FastCNMF programmes, evidence-based annotation, tumor CNV review, large-atlas scVI annotation and scPoli gating, spatial deconvolution, cell-cell communication, pseudobulk differential expression, abundance analysis, delivery, or maintenance of the SCOOP runtime and references.
---

# SCOOP

Use SCOOP as one skill with one managed CLI. Keep biological decisions in the
Skill and references; keep deterministic computation in the installed runtime.

## Runtime boundary

Resolve `SCOOP_HOME` from a non-empty absolute environment value. Otherwise use
`~/.scoop`. Require the directory to exist before analysis. Invoke only:

```bash
"${SCOOP_HOME:-$HOME/.scoop}/bin/scoop" ...
```

Do not use a repository script, host Python, Conda executable, or a bare
`fast*` command when the managed CLI is available. Run the relevant capability
check before computation:

```bash
"${SCOOP_HOME:-$HOME/.scoop}/bin/scoop" doctor --json
```

Read [runtime](references/runtime.md) when installing, diagnosing dependencies,
selecting an environment profile, or resolving paths.

## Operating contract

- Keep raw data read-only and retain raw counts in `layers["counts"]`.
- Keep project data and run outputs in the caller's project, not in SCOOP Home.
- Use semantic stage names. Do not create order-prefixed files or directories.
- Write large tables externally; keep only stable per-cell state and registry
  references in H5AD.
- Require every stage to emit a receipt with inputs, outputs, hashes, quality,
  software versions, warnings, and review state.
- Consume only committed upstream artifacts. Invalidate descendants when an
  input, configuration, label table, or approval-bound hash changes.
- Never silently switch algorithms. Record fallbacks in the decision log.
- Do not let an Agent directly mutate H5AD annotation fields. Export evidence,
  validate structured decisions, then use the deterministic committer.

Read [workflow contract](references/workflow-contract.md) for stage dependencies,
resume rules, artifacts, and hard gates. Read [design principles](references/docs/design-principles.md)
and [storage structure](references/docs/project-structure.md) when creating a run.

## Route the workflow

Choose the narrowest workflow that satisfies the request:

| Request | Route | Read before execution |
| --- | --- | --- |
| Validate input, QC, embedding, programmes, annotation, delivery | `standard` | [modules](references/docs/modules-and-outputs.md), [QC](references/stages/qc.md), [core](references/stages/core.md), [programmes](references/stages/programs.md), [annotation](references/stages/annotation.md) |
| Tumor annotation or malignant/normal separation | `tumor` | Standard references plus [tumor rules](references/markers/tumor-general.md) and the matching tissue reference |
| Multi-study or population atlas with scVI/scPoli | `atlas` | [large atlas](references/docs/large-atlas-integration.md), [annotation workflow](references/markers/annotation-workflow.md) |
| Spatial deconvolution or spatial graph | `spatial` | [spatial stage](references/stages/spatial.md), [parameter rules](references/docs/parameter-tuning.md) |
| Candidate cell-cell communication | `ccc` | [CCC stage](references/stages/ccc.md), [modules](references/docs/modules-and-outputs.md) |
| Marker, condition DE, abundance, or survival | `differential-expression` | [DE stage](references/stages/differential-expression.md), [modules](references/docs/modules-and-outputs.md) |

Create a plan before expensive computation:

```bash
"${SCOOP_HOME:-$HOME/.scoop}/bin/scoop" plan \
  --workflow <workflow> \
  --config <project>/runs/<run-id>/config/run.yaml
```

Execute only ready stages. Stop at any receipt with `review_required=true`.

## Annotation

Read [annotation workflow](references/markers/annotation-workflow.md) and
[naming strategy](references/markers/naming-strategy.md). Load only the relevant
marker references:

- [major cell types](references/markers/major-cell-types.md)
- [T cells](references/markers/t-cells.md)
- [myeloid](references/markers/myeloid.md)
- [fibroblast](references/markers/fibroblast.md)
- [tumor general](references/markers/tumor-general.md)
- [BRCA](references/markers/tumor-brca.md)
- [CRC](references/markers/tumor-crc.md)
- [GBM](references/markers/tumor-gbm.md)

Treat marker, negative-marker, programme, reference, CNV, neighbourhood, sample
distribution, and QC evidence as separate evidence classes. Never accept a
label from one method alone. Keep state programmes separate from cell identity.

## Atlas gate

Treat the scVI stage as representation learning followed by graph construction,
clustering, full annotation, confidence assignment, and a user-facing report.
scVI itself does not predict biological labels.

The report must show subtype distributions, graph-neighbouring groups,
positive/negative markers, reference conflicts, cluster-specific cNMF identity
and functional programmes, QC risks, confidence, and the exact cells proposed
for scPoli supervision.

Do not run scPoli until the managed CLI verifies an approved annotation gate.
The approval must bind the annotation report, label table, input H5AD fingerprint,
and labeled indices by SHA-256. Reject missing, rejected, incomplete, or stale
approval. Exclude legacy, scANVI, medium-confidence, low-confidence, ambiguous,
and suspicious cells from supervised indices.

## Delivery

Return the requested biological or analytical result in the final response;
do not return only internal paths. Also report the run ID, committed artifact,
review state, important warnings, and reproducibility manifest path.
