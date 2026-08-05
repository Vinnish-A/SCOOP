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

Optional branches do not change the dependencies of committed stages. Delivery
must select the newest committed H5AD from the actual dependency graph; it must
not guess by numbered filename or fall back from committed annotation to an
evidence-only artifact.

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

- QC removal above the configured review threshold pauses downstream work.
- Unstable embedding or clustering pauses annotation.
- Annotation decisions with unresolved conflicts pause commit.
- Malignant labels require tumor evidence or explicit review.
- scPoli requires a current approved annotation gate bound to exact artifacts.
- Delivery requires every mandatory ancestor to be committed.
