# SCOOP Fast & Annotation Decision Architecture

Target repository path: `docs/05_annotation_fast_architecture.md`

Status: draft architecture for implementation

## 0. Executive Summary

SCOOP should be structured as an evidence-native SOP rather than a loose toolkit.

The target architecture has four layers:

1. **SOP Workflow Layer**: run directories, module ordering, H5AD state, file registry, decision log.
2. **Fast Compute Layer**: deterministic engines such as FastCore, FastDE, FastCNMF, and FastCNVpy.
3. **Evidence & Skill Layer**: curated marker/state/tumor/naming rules, versioned as retrievable skills.
4. **Annotation Decision Layer**: schema-bound evidence fusion, validation, and commit into H5AD.

The Agent/AI layer must not directly mutate H5AD or invent labels. It may only propose structured annotation decisions that pass deterministic validators.

MCP is a future control-plane interface. Large matrices and long tables must stay in run artifacts; MCP tools should pass paths, manifests, and registry references, not matrix payloads.

---

## 1. Current Repository Grounding

Existing architecture already contains the right primitives:

- `README.md` defines SCOOP as an executable, auditable SOP rather than a tool collection.
- `docs/00_design_principles.md` says each module has one responsibility and annotation only fuses evidence.
- `docs/02_modules_and_outputs.md` defines the seven modules and already includes `04b_tumor_fastcnvpy`.
- `scripts/README.md` defines scripts as Agent execution interfaces.
- `src/scsp_agent_sop/decision_log.py` already supports decision evidence, fallback metadata, and human review flags.
- `scripts/04_annotation_markers.py` exports annotation evidence but intentionally does not assign final labels.
- `src/scsp_agent_sop/annotation.py` creates an evidence template with `review_required` defaults.
- `pyproject.toml` exposes `fastcore`, `fastde`, `fastcnmf`, and `fastcnvpy` as CLI entry points.

The missing layer is a deterministic **Annotation Decision Layer** that consumes evidence and skills, validates proposed labels, commits accepted labels, and logs the decision.

---

## 2. Target Layering

```text
SCOOP/
  scripts/
    04_annotation_markers.py          # existing evidence export
    04b_tumor_fastcnvpy.py            # existing tumor CNV evidence gate
    04c_annotation_decide.py          # new: produce/validate decision drafts
    04d_annotation_commit.py          # new: commit accepted decisions

  src/
    scoop_fast/                       # new: unified Fast engine contracts
      __init__.py
      artifact_bundle.py
      engine_spec.py
      quality_gate.py
      registry.py

    scsp_agent_sop/
      annotation_decision/            # new: evidence-native annotation layer
        __init__.py
        evidence_bundle.py
        decision_schema.py
        skill_models.py
        skill_retriever.py
        program_sanitizer.py
        validator.py
        committer.py

  configs/
    annotation_decision_schema.yaml   # new: allowed decision fields/gates
    default_run.yaml                  # update annotation/tumor config

  docs/
    05_annotation_fast_architecture.md
    skills/README.md                  # new: skill package convention

  tests/
    test_scoop_fast_contracts.py
    test_annotation_decision_schema.py
    test_annotation_validator.py
    test_program_sanitizer.py
```

---

## 3. Fast Compute Layer

### 3.1 Principle

Fast engines are deterministic compute units. They should not decide biological labels.

Each Fast engine should keep its independent package:

- `fastde`
- `fastcore`
- `fastcnmf`
- `fastcnvpy`

SCOOP should add a small contract layer under `src/scoop_fast/` that standardizes how SOP scripts and future MCP tools call these engines.

### 3.2 EngineSpec

`EngineSpec` is metadata about a Fast engine.

Required fields:

```python
@dataclass(frozen=True)
class EngineSpec:
    engine_id: str
    task_type: str
    version: str
    input_schema: str
    output_schema: str
    consumes: tuple[str, ...]
    produces: tuple[str, ...]
    default_cli: tuple[str, ...]
    writes_h5ad_fields: tuple[str, ...]
    writes_external_artifacts: tuple[str, ...]
    quality_gates: tuple[str, ...]
```

Initial engine IDs:

```text
fastde.markers
fastde.pseudobulk_deseq2
fastcnmf.programs
fastcnvpy.tumor_pooled
fastcore.preprocess
fastcore.core_pipeline
fastcore.graph_embed_cluster
fastcore.quality_compare
```

### 3.3 ArtifactBundle

All Fast engines should report results using the same lightweight bundle:

```python
@dataclass(frozen=True)
class ArtifactBundle:
    schema_version: str
    engine_id: str
    task_type: str
    run_id: str
    status: str
    inputs: dict[str, Any]
    outputs: dict[str, Any]
    quality: dict[str, Any]
    timings: dict[str, float]
    registry_patch: dict[str, Any]
    decision_log_patch: dict[str, Any]
    review_required: bool = False
    review_reason: str | None = None
```

This bundle is a control-plane summary. It must not contain full matrices or long result tables.

### 3.4 Registry

`src/scoop_fast/registry.py` should expose:

```python
list_engines() -> list[EngineSpec]
get_engine(engine_id: str) -> EngineSpec
```

Do not rewrite existing algorithms in this PR.

---

## 4. Annotation Decision Layer

### 4.1 State Machine

Annotation should become a stateful evidence workflow:

```text
RAW
  -> QC_DONE
  -> CORE_DONE
  -> CLUSTER_SELECTED
  -> MARKER_EVIDENCE_READY
  -> PROGRAM_EVIDENCE_READY
  -> BROAD_ANNOTATION_PROPOSED
  -> TUMOR_CNV_READY_OR_NOT_REQUIRED
  -> SUBCLUSTER_ANNOTATION_PROPOSED
  -> ANNOTATION_VALIDATED
  -> ANNOTATION_COMMITTED
  -> FINAL_MARKERS_READY
  -> ENRICHMENT_QC_READY
```

Invalid transitions must be rejected in future workflow code. The first implementation can encode these states as constants and validators without a full workflow engine.

### 4.2 Evidence Bundle

`AnnotationEvidenceBundle` is the input to AI/human annotation decisions.

Required fields:

```python
@dataclass(frozen=True)
class ClusterEvidence:
    cluster_id: str
    cluster_key: str
    n_cells: int
    n_samples: int | None
    marker_refs: tuple[str, ...]
    top_markers: tuple[str, ...]
    anti_marker_warnings: tuple[str, ...]
    nmf_refs: tuple[str, ...]
    dominant_programs: tuple[str, ...]
    program_warnings: tuple[str, ...]
    reference_refs: tuple[str, ...]
    cnv_refs: tuple[str, ...]
    cnv_summary: dict[str, Any]
    qc_warnings: tuple[str, ...]
```

```python
@dataclass(frozen=True)
class AnnotationEvidenceBundle:
    schema_version: str
    run_id: str
    organism: str
    tissue: str
    disease_context: str | None
    is_tumor: bool
    cluster_key: str
    clusters: tuple[ClusterEvidence, ...]
    file_registry: dict[str, Any]
```

### 4.3 Annotation Decision Schema

AI/human output must use a strict schema.

```python
@dataclass(frozen=True)
class AnnotationDecision:
    schema_version: str
    run_id: str
    cluster_id: str
    cluster_key: str
    parent_label: str
    canonical_label: str
    cell_state: str | None
    functional_modifier: str | None
    final_label: str
    confidence: str
    evidence_refs: dict[str, tuple[str, ...]]
    positive_markers: tuple[str, ...]
    negative_markers_absent: tuple[str, ...]
    conflicts: tuple[str, ...]
    review_required: bool
    reason: str
```

Allowed confidence values:

```text
high
medium
low
```

### 4.4 Validator Rules

`validate_annotation_decision(decision, evidence_bundle, skills)` must return a structured result, not raise for biological disagreement unless the schema is invalid.

Required validation checks:

1. The cluster exists in the evidence bundle.
2. `confidence` is one of `high`, `medium`, `low`.
3. `evidence_refs` references known artifacts or declared external evidence.
4. `final_label` must not be empty.
5. Program-dominated labels such as ribosomal, mitochondrial, stress, or cell-cycle should normally be written as `cell_state` or `functional_modifier`, not as a canonical cell type.
6. If `evidence_bundle.is_tumor` is true and the decision uses malignant/tumor/cancer labels, the decision must cite CNV evidence unless it is marked `review_required`.
7. If marker conflicts are present, confidence cannot be `high`.
8. If `n_cells` is below threshold or `n_samples` is too low, confidence cannot be `high`.
9. If the validator cannot verify enough evidence, it must force `review_required = True`.

### 4.5 Committer

`commit_annotation_decisions` is deterministic.

It should:

- load validated decisions;
- map clusters to cells through `adata.obs[cluster_key]`;
- write H5AD obs fields:
  - `cell_type_lvl1`
  - `cell_type_lvl2`
  - `cell_type_lvl3`
  - `cell_state`
  - `annotation_confidence`
  - `annotation_status`
- write external tables:
  - `annotation_decisions.tsv`
  - `annotation_validation.tsv`
  - `merge_split_log.tsv`
- register files through existing storage helpers;
- append decision logs.

The committer must not call an LLM.

---

## 5. ProgramSanitizer

Fine-grained annotation requires separating identity programs from nuisance/state programs.

`ProgramSanitizer` should classify NMF programs as:

```text
identity
ribosomal
cell_cycle
stress
mitochondrial
hemoglobin
ambient
mixed
unknown
```

Initial heuristic implementation:

- ribosomal: high fraction of genes starting with `RPL` or `RPS`;
- mitochondrial: high fraction of genes starting with `MT-`;
- hemoglobin: high fraction of `HBA`, `HBB`, `HBD`, `HBE`, `HBG`, `HBM`, `HBQ`, `HBZ`;
- cell_cycle: overlap with configured S/G2M gene sets or known genes such as `MKI67`, `TOP2A`, `UBE2C`, `HMGB2`, `PCNA`;
- stress: overlap with heat-shock/immediate early genes such as `HSPA1A`, `HSPA1B`, `JUN`, `FOS`, `DUSP1`;
- ambient/mixed: multiple incompatible lineage markers or high hemoglobin/mitochondrial signature outside relevant lineages.

Output table:

```text
program
program_class
action
top_genes
matched_genes
reason
```

Allowed actions:

```text
keep_for_identity
state_only
exclude_from_identity
review
```

---

## 6. Tumor CNV Gate

Tumor handling is mandatory after broad lineage annotation.

Rule:

```text
IF run.tissue/disease context indicates tumor:
  AFTER broad major lineage annotation:
    RUN FastCNVpy pooled tumor evidence
```

FastCNVpy evidence fields:

```text
fastcnv_reference_pool
fastcnv_cnv_fraction
fastcnv_normal_threshold
fastcnv_tumor_evidence
```

Decision validator rule:

- A malignant/tumor label must cite CNV evidence.
- If CNV is unavailable or indeterminate, the label must be downgraded to candidate/review, for example:
  - `Tumor-candidate epithelial cell`
  - `CNV-indeterminate malignant-like epithelial cell`

---

## 7. Skill Package Convention

A skill is a versioned biological rule package. It is not a free prompt.

Proposed layout:

```text
skills/
  human_gbm_tumor_v1/
    skill.yaml
    markers.tsv
    anti_markers.tsv
    states.yaml
    programs.yaml
    naming_rules.yaml
    ontology_map.tsv
    conflict_rules.yaml
    examples/
      cluster_evidence_001.json
      accepted_decision_001.json
```

Initial code should only define the schema and loader. Do not build a large database in this PR.

### 7.1 Skill YAML

```yaml
skill_id: human_gbm_tumor_v1
version: 1.0.0
species: human
tissue: brain
disease_context: glioblastoma
lineage_scope: tumor_and_microenvironment
gene_symbol_namespace: HGNC
sources: []
naming_policy:
  uncertain_format: "{phenotype} - {function}"
  use_cell_state_for:
    - cycling
    - ribosomal_high
    - stress_high
```

---

## 8. MCP Direction

Do not implement MCP in the first PR.

The architecture must stay MCP-ready:

- tools receive paths/config/run IDs;
- tools return `ArtifactBundle`;
- no matrix payloads in tool responses;
- no arbitrary shell execution;
- long tables remain external artifacts and are referenced by registry entries.

Future namespaces:

```text
scoop-compute-mcp
scoop-knowledge-mcp
```

---

## 9. Implementation Plan

### Phase 0: Documentation and Config

Add:

- `docs/05_annotation_fast_architecture.md`
- `docs/skills/README.md`
- `configs/annotation_decision_schema.yaml`

Update:

- `configs/default_run.yaml`
  - add `annotation.decision`
  - add `annotation.program_sanitizer`
  - add `tumor_fastcnv`
- `scripts/README.md`
  - include `04b_tumor_fastcnvpy.py`, `04c_annotation_decide.py`, and `04d_annotation_commit.py` in the intended flow.

### Phase 1: Fast Contract Skeleton

Add `src/scoop_fast/` with:

- `EngineSpec`
- `ArtifactBundle`
- simple registry for current Fast engines
- quality gate helpers

Tests:

- registry returns expected engine IDs;
- bundles round-trip through JSON.

### Phase 2: Evidence Bundle Builder

Add `src/scsp_agent_sop/annotation_decision/evidence_bundle.py`.

Implement:

- build bundle from AnnData obs, file registry, cluster key, marker table path if present;
- include CNV summaries when FastCNVpy fields exist;
- include NMF summary when `dominant_nmf_program` or `X_nmf_usage` exists;
- degrade gracefully when optional evidence is absent.

Tests:

- evidence bundle builds on tiny AnnData;
- missing optional fields do not crash.

### Phase 3: ProgramSanitizer

Add `program_sanitizer.py`.

Implement lightweight gene-prefix/marker-overlap heuristics.

Tests:

- ribosomal top genes classify as `ribosomal`;
- cell-cycle top genes classify as `cell_cycle`;
- mixed/unknown programs return review or unknown.

### Phase 4: Decision Schema and Validator

Add:

- `decision_schema.py`
- `validator.py`

Implement:

- JSON load/dump helpers;
- dataclass validation;
- tumor CNV gate;
- confidence downgrade rules;
- forced review rules.

Tests:

- malignant tumor label without CNV evidence fails or forces review;
- marker conflict prevents high confidence;
- invalid confidence is rejected.

### Phase 5: Committer and Scripts

Add:

- `committer.py`
- `scripts/04c_annotation_decide.py`
- `scripts/04d_annotation_commit.py`

Behavior:

- `04c` builds evidence bundle and writes a decision template JSON/TSV. It must not call an LLM.
- `04d` reads user/AI-edited decisions, validates them, commits accepted fields into H5AD, writes audit tables, and logs decisions.

Tests:

- committing decisions writes expected obs fields;
- invalid decisions are not committed as accepted;
- output files are registered.

### Phase 6: Tumor Flow Integration

Update documentation and config so tumor samples run:

```text
04_annotation_markers.py
04b_tumor_fastcnvpy.py
04c_annotation_decide.py
04d_annotation_commit.py
```

Do not reimplement FastCNVpy.

### Phase 7: Final Marker Recalculation Hook

Add documentation or lightweight stub for recomputing markers after final annotation.

Do not implement full enrichment in the first PR unless trivial.

### Phase 8: MCP Adapter

Future PR only.

---

## 10. Non-goals for the First Implementation PR

Do not:

- rewrite FastCNMF internals;
- rewrite FastCNVpy;
- implement a live LLM annotation agent;
- implement MCP server;
- introduce a large marker database;
- store long marker/NMF/CNV/DE tables inside H5AD;
- allow arbitrary shell commands through the new layer;
- change existing default biological outputs unless a test covers the change.

---

## 11. Acceptance Criteria

The first architecture implementation PR is accepted if:

1. `pytest -q` passes.
2. New unit tests cover Fast contracts, annotation decision schema, validator, and committer.
3. Existing scripts remain backward compatible.
4. The new `04c` script can generate a decision template from a small H5AD.
5. The new `04d` script can commit a valid decision file into H5AD.
6. Tumor label validation requires CNV evidence or review.
7. Large result tables remain external artifacts.
8. Every commit operation writes decision-log records.

---

## 12. Codex Implementation Scope Recommendation

Implement only Phases 0-5 in the first PR.

Phase 6 should be mostly documentation/config wiring unless tests are easy.
Phase 7 should be a documented hook, not a full enrichment implementation.
Phase 8 should remain deferred.
