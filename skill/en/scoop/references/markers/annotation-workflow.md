# General Single-Cell Annotation Workflow

When `identity_scope` is not declared, record the accepted granularity as `broad` only. Fine-level adjudication must explicitly declare
`identity_scope: fine`; `accepted` or “reviewed RNA phenotype” alone does not establish that fine subtypes have been resolved.
When migrating old tables, check the evidence actually reviewed; do not infer fine-level resolution from names, confidence, or a single status field.

Version: v1.0-refined  
Goal: complete scRNA-seq cell annotation quickly, robustly, and simply, preserving evidence from every step in reviewable tables.

---

## 1. Minimum inputs

```yaml
inputs:
  adata: normalized scRNA-seq object
  clusters: reviewed clustering outputs and within-parent neighbourhoods
  cluster_deg: differential genes for each cluster
  cnmf_programs: batch-handled parent cNMF / FastCNMF programs and usage, required for formal fine review, not broad review
  cnv_result: FastCNVpy result for tumor-relevant compartments; supporting evidence, never required
  project_context:
    species: human | mouse | other
    tissue: <tissue>
    disease: tumor | normal | inflamed | other
    cancer_type: BRCA | CRC | GBM | other
    sample_info: patient, region, treatment, pathology if available
```

---

## 2. Overall workflow

For an integrated atlas, annotation must center on unsupervised clustering and actual gene-level
DEG within the current integration scope; old studies and per-source identities are references only. Naming each source and concatenating them is not post-integration reannotation,
and old-label allowlists must not exclude competing lineages missed by the original study. First verify reusable
QC/clustering/CNV/programme checkpoints individually, then fill gaps; see [multimodal preparation and reuse](../stages/multimodal.md).
Once clustering, its inputs, and measurement scope have been verified as reliable, directly compute conventional cluster,
subcluster, and merged final-label gene-level DEG on the actual expression layer; source/batch stratification is an optional diagnostic
for specific questions, not an annotation prerequisite. Check source composition, denominators, and interpretive limits from batch confounding;
cell-level differences do not replace independent replication at the patient/condition level.
For key markers absent from HVGs/shared columns, return to the source’s measured gene axis for verification; a missing column is not a negative result.
Marker means, detection rates, and specificity rankings are descriptive evidence and must not be presented as DEG with actual P/q values.
Native `annotate evidence` exports statistical tests over the full input gene axis; the COSG option separately saves an auxiliary
specificity table and does not replace DEG. Review templates inherit registered cluster_markers references; identity decisions without
references to clustering differential evidence remain review_required. The commit interface verifies registered reference files
and hashes, while the analyst still adjudicates their biological meaning, denominators, and comparison scope.

Preserve original signed test statistics; for example, calculate Wilcoxon/MWU Z from U and the actual tie-adjusted variance,
rather than back-calculating Z from clipped P after P underflows to 0. The latter creates many ties among strongly differential genes, making top
lists depend on gene order. Verify P/q, U, effect sizes, detection rates, and ranking statistics separately; for export
problems, derive missing statistics from verified U/expression axes without refitting valid integration/graphs/clustering.

First-level broad classes are reviewed, biologically meaningful cell classes such as T/NK, rather than numeric clusters
or unreviewed mergers of numeric groups. The main process may adjudicate broad classes using reliable clustering, conventional DEG, and known tissue
positive/negative markers; broad classes do not require NMF as a prerequisite. Numeric groups may be merged or split according to evidence,
with boundaries reviewed again and competing and evidence-pending cells retained; proving every cell’s purity is not required.

Before formally starting fine-subgroup analysis, have a reviewed broad biological parent class and usable metadata, batch, QC,
and measurement contracts, and complete that parent class’s programme NMF with appropriate actual expression-level or within-model
batch handling. Once the parent is reviewed, NMF may run in parallel with preparatory within-class dimensionality reduction, clustering, and conventional DEG;
while gaps remain, this work serves only to supplement evidence or explore candidates, and must not be called completed formal fine-subgroup analysis.
Before formal fine-subgroup adjudication, all evidence in section 5 must also be collected; if any item is missing, record
`evidence_pending`, actively obtain evidence, and do not submit candidates as final identities. Retain UIDs, sources, and competing context for the full declared parent scope; maximum programme
usage is not an identity classifier. Exploratory NMF on a numeric scope does not automatically satisfy the biological parent-class requirement.

When subdivision stability is insufficient, retain results for every seed and resolution, distinguishing stable broad lineages, uncertain boundaries,
and continuous states; do not lower declared thresholds or force names merely to deliver. First check actual measurements and representation
features: when key measured genes were not included in atlas-wide HVGs, a comparable feature space may be refitted
on the newly integrated whole-parent scope; this differs from unmeasured genes on the native axis, and neither gap can be solved by zero filling or old-label filtering.

When subdividing an inherited atlas-wide representation, first quantify how well actual fine-type gene signals are retained within the scope, and assess this together with
actual DEG, supplementary evidence for missing native columns, and sample structure. Well-preserved representations can reconstruct local neighborhoods
for the whole new parent population; refit the relevant scope only when there is specific information loss, and do not discard the entire
representation by default because of one low ARI. Record insufficient optimization under a fixed iteration budget, continuous states, and unreliable identity boundaries separately.

Adjudicate labels for identical UIDs across scopes using actual expression/DEG and competing evidence jointly; do not overwrite them directly based on scope names,
old routing, or maximum confidence. A reviewed recursive subgroup replaces only the parent
scope it actually refines, retaining its ancestor models/row numbers and all UIDs; conflicts in other independent scopes still require review. Final
cell tables and merged-label DEG are the actual deliverables; naming notes or a mapping proposal alone do not constitute a completed commit.

```text
Step 0. QC / ambient RNA / doublet review
Step 1. Select final clustering resolution
Step 2. Compute DEG for every cluster
Step 3. Main agent reviews actual clustering, conventional DEG, tissue markers and neighbourhoods to accept meaningful biological broad classes; NMF is not required
Step 4. Freeze the accepted parent UID scope, technical contracts and mixed, conflicting or evidence-pending boundaries
Step 5. Prepare parent batch-handled NMF and within-parent embedding/clustering/DEG, optionally in parallel
Step 6. Start formal fine-subtype analysis only after the parent, technical prerequisites and parent NMF are complete; adjudicate only after all section 5 evidence is complete
Step 7. Review tumor-specific CNV/SNV and context separately from lineage/state when relevant
Step 8. Attach programmes as functional overlays; programme argmax is not an identity classifier
Step 9. Merge redundant reviewed subclusters into canonical final labels when appropriate
Step 10. Recompute DEG per final annotated subtype
Step 11. Export marker tables, annotation evidence table, and quick enrichment report
```

Above, cNMF/FastCNMF denotes the programme task; the actual execution engine is determined by [programmes](../stages/programs.md)
and the receipt. Do not present native `sklearn_multistart_nmf` as standalone `fastcnmf`.

---

## 3. Decision branches

### 3.1 Non-tumor samples

```text
credible clustering + conventional DEG + tissue markers -> reviewed biological broad class
  -> batch-handled parent NMF and preparatory within-parent clustering/DEG (may be parallel)
  -> completed parent/technical/NMF prerequisites -> formal subtype analysis
  -> all section 5 evidence + main-agent adjudication -> reviewed fine label
```

### 3.2 Tumor samples: review malignancy evidence independently

```text
broad annotation
  -> identify epithelial / parenchymal / glial / other candidate malignant compartments
  -> CNV/SNV, normal reference and patient/tissue-context review
  -> declared malignancy_status (lineage markers + DEG decide; CNV/SNV support) separate from differentiation/state
  -> cancer-type references support candidate review; malignant subtype names require supported malignancy
  -> normal tissue skill for normal-reference/context-supported parenchymal cells; CNV-neutral alone does not route cells to normal
```

#### Case A: normal parenchymal and tumor cells are already separated

FastCNVpy provides supporting evidence:

```text
tumor-like programme + concordant CNV/context -> reviewed malignant support
coherent normal lineage + matched normal-reference + tissue-context support -> normal; CNV neutrality alone is insufficient
```

#### Case B: normal parenchymal and tumor cells are not separated

Distinguish origins using joint evidence:

```text
same marker-defined parenchymal/epithelial cluster
  -> CNV abnormal cells: supporting evidence; review reference, programme and context
  -> CNV neutral cells: inspect normal atlas, tumor programmes, patient structure and tissue context; missing evidence is pending, while reviewed conflicting evidence may be unresolved
  -> CNV weak or absent with concordant non-CNV tumor evidence: retain probable_malignant as an option
  -> unresolved multi-source conflict: ambiguous; record alternatives and confidence
```

---

## 4. Broad annotation rules

Broad groups should be biologically meaningful, evidence-supported first-level cell classes (such as T/NK); numeric groups may be
merged or split according to actual evidence. The table below gives examples of lineage clues; a single marker or numeric group does not automatically become
a final class. Actual broad-class granularity must be determined by tissue context and competing evidence:

```yaml
broad_type:
  immune: PTPRC-positive compartments
  epithelial: EPCAM/KRT/CDH1-positive compartments
  stromal: COL1A1/DCN/LUM/PDGFRA-positive compartments
  endothelial: PECAM1/VWF/KDR/CLDN5-positive compartments
  mural: RGS5/CSPG4/PDGFRB/MYH11-positive compartments
  erythroid_or_platelet: HBB/HBA/PPBP/PF4-positive compartments
# Missing evidence is tracked in workflow progress, not as a broad_type.
```

Do not assign excessively fine subtype names at the broad-group stage, for example:

```text
Do not directly call T_CD8_TEXterm during initial clustering of all cells.
Do not directly call Mac_TAM_SPP1_TREM2 before extracting and reclustering myeloid cells.
Do not directly call MaligEpi_LuminalB in a tumor epithelial group without a reviewed malignancy decision.
```

---

## 5. Fine annotation rules

Fine subgroups are biologically meaningful populations, rather than numeric IDs. Before formal analysis, have a reviewed biological parent class,
metadata/batch/QC/measurement contracts, and batch-handled parent NMF. Within-class
dimensionality reduction, clustering, UMAP neighborhoods, and conventional DEG for candidate partitions may be generated in parallel with NMF as preparatory evidence.
Before formal adjudication, the main process must personally read the actual results and collect:

1. The reviewed broad parent class and its UID scope, actual metadata, batch, QC, measurement comparability, and unknowns;
2. Parent NMF spectra, usage, stability, actual batch handling, and input scale;
3. Within-class dimensionality reduction, clustering results, and UMAP distances/neighborhoods, assessing boundaries together with actual expression; two-dimensional UMAP distances alone cannot determine lineage;
4. Conventional subcluster DEG, positive/negative and competing-lineage markers, tissue knowledge, and available relevant context;
5. Cross-sample support, denominators, and checks for batch/condition confounding, conflicts, or missing measurements; not every subgroup must occur in every experimental condition, since condition specificity may be real biology;
6. Separately review CNV/SNV or tumor references when tumor-relevant; use an appropriate reference atlas as supporting evidence when available.

If any required item is missing, mark workflow status `evidence_pending` and actively retrieve original sources, tissue
markers, and verifiable literature, or run necessary supplementary analyses; record which evidence is unavailable, unmeasurable, or temporarily
inexecutable and why, and continue feasible work. Candidates and gaps must be explicitly recorded; do not formally name or submit them as
complete. Use identity `unresolved` only when sufficient evidence has been obtained and the main process’s review still cannot adjudicate competing, mixed, or conflicting evidence;
do not use unresolved to conceal evidence that has not been retrieved or completed.
`evidence_pending` is a work-progress or evidence-gap-table status; do not force a new enum into the existing annotation
schema or present it as a cell identity. Cell-level marker DEG describes identity and does not replace independent replication and
corresponding statistical inference for patient/condition effects; conventional cluster DEG need not be forcibly stratified by source or batch.

Uses of cNMF / FastCNMF:

```text
identity-disambiguation: after removing uninformative dominant axes such as cycling/ribosomal/stress, inspect the remaining core identity programmes
functional-overlay: attach IFN, hypoxia, cytotoxic, APC, remodeling, and similar signals as state suffixes
cluster-merge: if two small groups differ only in cNMF functional strength but share the same historical phenotype, they may be merged under one name
rare-state-detection: give small numbers of special states low/medium-confidence state labels, rather than creating new “cell types”
```

---

## 6. Subagent annotation rules

**The main process schedules and reviews; subprocesses may run analyses, retrieve and integrate evidence, and propose candidates, but final annotation must be adjudicated by the main process personally reading actual DEG, clustering, UMAP distances/neighborhoods, and tissue markers; for fine subgroups, it must also personally read actual parent NMF and all evidence in section 5. Do not merely accept subagent labels. Do not rely on mechanical classification or automatic label propagation with SingleR, CellTypist, scANVI, or similar tools to determine final names; atlas matching is supporting material only. Record provenance, alternative identities, reasoning, and confidence for each group, and commit after structured validation. Treatment response/survival outcomes do not participate in identity naming. CNV neutrality alone also does not establish normality; normal identity and tissue-context support are still required.**

Subprocesses may organize candidate evidence in the following order; the main process makes the final judgment after independently reading the actual evidence:

```text
1. Check broad lineage gate.
2. Reject incompatible lineages using exclusion markers.
3. If tumor-relevant, review malignancy_status on lineage markers and DEG, with FastCNVpy/SNV/reference as support.
4. Search matching historical phenotypes in the relevant marker skill.
5. If a canonical historical phenotype is supported, use it.
6. If phenotype is uncertain but real cell identity is clear, use <real phenotype + function>.
7. If two small clusters map to the same real phenotype and only differ by generic cNMF programs, merge naming.
8. Record confidence and do_not_call labels.
```

---

## 7. Confidence

```yaml
high:
  required:
    - broad lineage marker set consistent
    - at least two to three core phenotype markers consistent
    - key exclusion markers absent or low
    - tissue / cancer / CNV / reference context consistent
medium:
  required:
    - broad lineage reliable
    - phenotype marker set partially consistent
    - some uncertainty in exclusion, CNV, or reference mapping
low:
  triggered_by:
    - single-marker call
    - ambient RNA or doublet risk
    - weak or conflicting malignancy evidence in a tumour-relevant lineage
    - conflicting historical phenotypes
    - very small cluster or low UMI
```

---

## 8. Final exports

After each annotation round, export:

```text
cluster_markers.tsv/parquet
subcluster_markers.tsv/parquet
final_label_markers.tsv/parquet
annotation_evidence.tsv
cnmf_program_summary.tsv
fastcnvpy_summary.tsv, if tumor
quick_enrichment_by_final_label.tsv
```

Recompute `final_label_markers` after final naming and merging; do not simply reuse original cluster DEG.


### Non-empty validation of exported review tables

When transferring DE statistics into per-group review tables, compare directly against the actual tested gene axis and original
DEG rows. Explicitly distinguish unmeasured genes, measured but undetected genes, genes untested because of insufficient denominators, and genes absent from
the requested panel; these are not interchangeable. For gene-indexed tables, check the gene index;
do not mistakenly inspect DataFrame column names and record actual measured markers as missing.

Validation must report the number of group×gene rows actually checked and the number of non-empty fields, with independent spot checks of
U, signed Z, P, BH q, mean FC, and detection rates for positive, negative, highly significant, and untested
records. An all/pass result on an empty set is not evidence that the review table was validated.
If an export-logic error hides original DE evidence from review, archive the invalid review version and
review all affected identities and states again; when original expression and tests are unchanged, do not rerun models just to fix a report.

### Review coverage and identity granularity

Structured decisions use `identity_scope: broad | fine | fine_unresolved | unresolved` (default `broad`; see the
[structured decision contract](#structured-decision-contract)). Keep missing evidence in work-progress/gap tables first and actively obtain it; do not submit unresolved merely because evidence is missing.
That status is only for adequate evidence that remains mixed or conflicting after review. Committing writes `annotation_reviewed`, `broad_annotation_status`, and
`fine_annotation_status` separately; reviewing only a broad lineage does not automatically review fine subgroups or confer normality.
A rejected new decision must revoke the corresponding cells’ old accepted eligibility, retaining old text for provenance.
Explicitly distinguish all records reviewed, broad lineage resolved, fine subgroups resolved, QC measured, and reference eligibility.

Evidence from cluster sampling may conceal small heterogeneous residual populations. Separately inspect competing lineages
and sampling coverage in the final unresolved UID subset; return to actual native features, unsupervised clustering, and all-gene DEG when systematic competing signals exist.
Do not describe missed checks, unrequested panels, or lineage genes absent from HVGs as intrinsic unidentifiability of the data.

### Structured decision contract

- **Evidence binding.** `annotate review` writes `annotation_evidence_bundle.json` (with `partition_sha256` over
  obs names and cluster-key values, and `deg_tables {ref: label_key, groups, sha256}`) and a template whose
  decisions carry `bundle_sha256`; the previous pair is archived byte-for-byte under `history/`. A DEG citation
  `ref` counts only when that table's `label_key` is the decision's cluster key and it tests the unit's group; cite
  another partition's group as `ref::group` (in the unit's own partition `::group` may name only the unit itself). A
  table registered without `label_key` binds no citation; re-register it with the partition it tests.
  `manual:`/`external:` refs remain declared, unbound DEG. Commit finds the bundle the decisions cite by its sha256
  (next to the decisions, in `annotation.decision.output_dir`, or their `history/`), and refuses a changed
  partition, a unit whose cell count differs from the bundle, and decisions whose bundle cannot be found. A legacy
  bundle without a partition hash is rebound from the input (`limited`); with no review bundle at all, commit builds
  one and reports the partition binding `not_assessed`.
- **Support.** `n_samples` counts samples with at least `annotation.decision.min_cells_per_sample` (default 20) cells
  in the unit (`n_samples_raw` keeps all). Composition keys come from `annotation.decision.{sample,donor,study}_key`,
  else `run.obs_keys`. `min_cells_high_confidence` and `min_samples_high_confidence` reach the validator.
- **Scope.** `fine_unresolved` accepts the broad lineage and records the fine level as reviewed but unresolved
  (`fine_annotation_status: unresolved`, no forced review). `unresolved` writes the labels but forces review.
- **Technical outcome.** `disposition: technical_exclusion` with `technical_kind: doublet | multiplet | ambient |
  low_quality | artefact` commits `annotation_status: excluded_technical`, broad/fine status `not_applicable`; the
  cells count as reviewed (`technical_excluded`) and remain deliverable. It needs unit DEG or a cited artefact
  registered under a QC/doublet schema (`qc*`, `scrublet*`, `doublet*`, `ambient*`, `core_quality*`); the sample-support
  rule for high confidence does not apply to an artefact.
- **Annotations without gating.** `notes` never gates; keep `conflicts` for genuine marker conflicts (they block
  high confidence). `parent_unit`/`unit_rule` record unit lineage in `merge_split_log.tsv`, which logs every
  decision as `commit | record_unresolved | technical_exclusion | review_required`.
- **Per-cell layer.** `annotate commit --cell-layer <table>` (`obs_name` plus any of `cell_confidence`, `cell_state`,
  `identity_source`) is bound as an input and applies only to cells whose unit labels were written. Per-cell
  confidence can only lower the unit confidence; `identity_source` fills `annotation_identity_source`.
- **Re-commit.** A unit re-decided without writing labels keeps its label text for provenance, but its
  `malignancy_status` returns to `not_assessed` and basis, technical kind and identity source are cleared.
- **Malignancy** is declared, never read from label words; see [tumor rules](tumor-general.md#4-malignancy_status-labels).
