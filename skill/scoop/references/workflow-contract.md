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
status: running | executed | failed | awaiting_review | prepared | skipped
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

`executed` means the command returned and its declared artifacts were recorded.
Read `execution`, `validation`, `interpretation` and measured `coverage`
independently. Legacy `committed` receipts imply none of the missing assessments.
Old receipts are archived under `state/stages/history/` before another attempt;
directory inventories do not become outputs of the new attempt.

Use `scoop status --run-root RUN --verify` to inspect declared artifact hashes,
missing/stale results and unassessed work. Hash integrity verifies files, not
biological correctness. There is no universal green “all analysis complete”.
Bind relevant upstream decisions/receipts as inputs when their replacement must
invalidate a dependent result; undeclared dependencies cannot be detected.

Project adapters can use the same result contract:

```json
{
  "stage": "spatial",
  "inputs": ["counts.h5ad", "reference_profiles.npz", "review.json"],
  "outputs": ["weights.parquet", "checks.json"],
  "execution": {"status": "completed"},
  "validation": {"status": "limited", "scope": "numerical_and_axis_checks"},
  "interpretation": {"status": "not_assessed"},
  "coverage": [{"scope": "fitted_spots", "unit": "observation_views", "total": 100, "processed": 95}],
  "limitations": ["Related fine components remain unresolved"]
}
```

Register with `scoop register-result --run-root RUN --manifest result.json`.
Relative artifact paths resolve against the manifest. Completed registrations
require actual outputs; `validation.status=passed` requires bound inputs and a
nonempty list of passed `checks`. These are explicit analyst declarations, not
independent adjudication by the registrar.

## Gates

- Review QC removal and unstable embeddings before dependent interpretation.
- Review annotation conflicts before commit; an explicit unresolved identity is
  a valid reviewed result, not a reason to stop every branch.
- Malignant labels require joint biological evidence/review; CNV neutrality does
  not establish normality.
- scPoli requires a current approved annotation gate bound to exact artifacts.
- Reviewed annotation delivery requires per-observation review coverage and a
  committed input. Explicitly unresolved reviewed cells remain deliverable.
  `scoop deliver --scope current-state` exports an earlier checkpoint deliberately;
  the default cannot silently fall back to QC/programme/evidence files.

Review pauses affected decisions/dependents, not independent work. Tool bypass
cannot replace missing biological evidence or supervision approval.
