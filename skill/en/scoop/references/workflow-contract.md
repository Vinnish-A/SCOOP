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
A result with `execution: not_run` (for example a tumour-CNV plan that cannot run)
records the stage `skipped` and keeps its `execution.reason`; a registered
`not_run` result is `skipped` too.
Read `execution`, `validation`, `interpretation` and measured `coverage`
independently. Legacy `committed` receipts imply none of the missing assessments.
Old receipts, including `atlas-approval`, are archived under
`state/stages/history/` before another attempt; directory inventories do not
become outputs of the new attempt. Each `spatial.rctd.output_label` keeps its own
receipt, `state/stages/spatial-<label>.json`.

The configuration input records the whole `run.yaml` sha256 (`config_hash`) as
provenance, plus the sections the stage reads (`config_sections`, for example
`run`/`keys`/`qc` for qc) and their hash. Staleness follows the section hash, so
editing an unrelated section does not mark the stage stale. Every stage also binds
the data file it read: `--input`, or its default input when the command declares
none.

`--output-format obs-parquet` (qc, tumor-cnv, annotate commit) writes the
stage's declared obs columns keyed by `obs_name`, whether or not a value changed,
and a `<stem>.json` receipt binding the table and its base input (path, sha256,
obs_name order). Annotate review and commit take their default input from what
the last `tumor-cnv` receipt declared, its sidecar or its H5AD, never from which
files exist; without a tumor-cnv receipt they read
`adata_annotation_evidence.h5ad`. A sidecar input, by default or `--input`, is
overlaid onto its verified base; a changed base or row order is refused.

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
- Malignancy is declared per decision (`malignancy_status`/`malignancy_basis`);
  lineage markers plus a valid DEG citation decide it. CNV is optional support and,
  when cited, must be a registered CNV artefact with assessed cells. CNV
  neutrality does not establish normality. See [tumor rules](markers/tumor-general.md#4-malignancy_status-labels).
- Annotation commit binds decisions to the reviewed content: the bundle's partition
  hash, unit cell counts, `bundle_sha256` and DEG tables that test the unit; see the
  [structured decision contract](markers/annotation-workflow.md#structured-decision-contract).
- scPoli requires a current approved annotation gate bound to exact artifacts.
- Reviewed annotation delivery requires per-observation review coverage and a
  committed input. Explicitly unresolved and `excluded_technical` reviewed cells
  remain deliverable.
  `scoop deliver --scope current-state` exports an earlier checkpoint deliberately;
  the default cannot silently fall back to QC/programme/evidence files. The atlas
  route writes no committed annotation, so deliver it with `--scope current-state`
  or commit reviewed decisions with `scoop annotate commit` first. The delivery
  schema alone names the retained obs/var/uns keys.

Review pauses affected decisions/dependents, not independent work. Tool bypass
cannot replace missing biological evidence or supervision approval.
