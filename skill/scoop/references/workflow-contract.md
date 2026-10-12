# Workflow contract

## Semantic workflows

```text
standard:
  validate -> qc -> core -> programs -> annotation-evidence
  -> annotation-review -> annotation-commit -> deliver

tumor:
  standard through annotation-evidence -> tumor-cnv
  -> annotation-review -> annotation-commit -> deliver

atlas:
  validate -> qc -> atlas-integrate -> programs -> atlas-annotate
  -> atlas-report -> atlas-approval -> atlas-scpoli
  -> atlas-validate -> deliver

spatial:
  annotation-commit -> spatial

ccc:
  annotation-commit -> ccc

differential-expression:
  annotation-commit -> de-prepare -> de-run
```

Use the narrowest route; these managed dependencies do not require every module.
Unusable interfaces permit validated `project-adapter` artifacts with reason,
method/environment, inputs/config, outputs/hashes and checks in existing records.
Do not fabricate managed approval/completion. Follow actual validated dependencies,
not numbered filenames; reviewed annotation cannot become evidence-only output.
Invalidate only descendants affected by changed inputs, labels or configuration.

## Stage receipt

Every stage writes an atomic JSON receipt containing:

```yaml
schema_version:
run_id:
stage:
status: planned | running | committed | failed | awaiting_review | invalidated
inputs:
input_hashes:
config_hash:
code_version:
outputs:
output_hashes:
quality:
warnings:
review_required:
review_reason:
software_versions:
started_at:
finished_at:
elapsed_seconds:
```

A committed receipt is reusable only when input, configuration, code and output
hashes still match. Failed work remains inspectable. Never erase it to simulate
a clean retry.

## Gates

- Review QC removal and unstable embeddings before dependent interpretation.
- Review annotation conflicts before commit; an explicit unresolved identity is
  a valid reviewed result, not a reason to stop every branch.
- Malignant labels require joint biological evidence/review; CNV neutrality does
  not establish normality.
- scPoli requires a current approved annotation gate bound to exact artifacts.
- Delivery requires complete, validated ancestors, managed or explicitly adapted.

Review pauses affected decisions/dependents, not independent work. Tool bypass
cannot replace missing biological evidence or supervision approval.
