# General Skill for Malignant / Normal Adjudication in Tumor Tissue

Version: v1.0-refined  
Scope: tumor samples in single-cell transcriptomics.  
Purpose: determine whether cells are malignant before entering a cancer-specific tumor skill.

---

## 1. General principle

```text
Establish malignancy first, then subtype within the cancer type.
```

In tumor tissue, normal parenchymal and tumor cells often share markers. For example:

- Both normal breast luminal cells and luminal breast cancer can express `EPCAM/KRT8/KRT18/KRT19/GATA3`.
- In CRC, normal crypt/absorptive/goblet cells and tumor epithelium can all express intestinal epithelial markers.
- In GBM, markers of normal astrocytes/OPCs/neurons overlap extensively with AC-like/OPC-like/NPC-like tumor markers.

Therefore, markers can establish only a “candidate origin,” not malignancy on their own.

### 1.1 Lineage markers decide; RNA-CNV supports

Identity is decided by the coherent lineage marker programme. Malignancy is assessed only among cells of a
tumour-candidate lineage, and it rests on convergent evidence: tumour programmes, aberrant co-expression,
expressed amplicon genes, patient-specific structure and reference mapping. Expression-inferred CNV is
supporting evidence in that review and is never the deciding criterion on its own. When lineage markers
and RNA-CNV disagree, the markers win and the discordance is recorded next to the label.

RNA-CNV is lineage- and reference-confounded. A lineage that is absent from the CNV reference, or whose
programme clusters on one chromosome (for example ECM genes on chr7), can produce apparent arm-level calls.
Summaries that depend on the reference scope (q99 exceedance, +7/−10 flags) move with the reference
composition. Before using CNV to separate cells, check it on externally validated cells of the same lineage.

A cell that passes a non-tumour lineage gate (endothelial, mural, fibroblast, immune) without co-expressing
tumour-lineage markers keeps that lineage whatever its CNV. Co-expression of a full non-tumour programme
with tumour-lineage markers is resolved per cell, as doublet/ambient or as an explicit `<A>_or_<B>_lowConfidence`
label. It is never resolved as transdifferentiation.

---

## 2. Where FastCNVpy is used

### 2.1 Broad-group stage

```text
broad cell types -> candidate compartment -> CNV/SNV/reference/context review -> malignant_status
```

Candidate malignant compartments typically include:

```text
epithelial/parenchymal/glial/neural crest-derived tumor cells
```

Immune, endothelial, fibroblast, and pericyte populations may be candidate normal references, but their complete identities, mixed origins, technical quality, and tissue context must be reviewed; old labels or CNV-low status do not directly confer reference eligibility.

### 2.4 Actual references, measurement scales, and checkpoints

Record both requested and effective reference UIDs/counts, both denominators for each sample,
minimum-reference parameters, excluded references and reasons, and the sample scope actually involved in gene selection,
baseline scaling, and thresholds. Passing an allowlist to the outer interface does not mean the backend
uses it in full; samples with few references may be excluded by the default minimum, potentially even changing gene-selection scope.
Inspect actual backend selection/return structures and retain unused references’ original positions in the complete query.
State whether the baseline is weighted by cell or by sample; the total number of references passed in cannot replace the actual denominator.

Bind checkpoint identity to count content/UID/gene order, sample key, reference
labels and actual UID hashes, minimum, samples actually used for gene selection, gene metadata/build/arm
policies, backend version/source code, and actual parameters. Effective parameters outside the configuration object must also be bound.
Create a new receipt when inputs or baselines change; retain old continuous results, but do not present them as the same fit.

Even after successful per-UID bridging, distinguish old/new count content and the shared-gene comparison scope. Continuous profiles from an old input version
may provide support with that version stated; they must not be relabeled as recomputation on current native counts or copied to uncovered cells
in the same group. Candidate and reference roles may overlap; candidate=False cannot
replace actual reference membership lists. Check how many patient×mode combinations new QC affects against the actual reference set difference,
rather than using it to invalidate all other valid fits. Record the review status of existing CNV evidence separately from final malignancy
assignment: retain undecided cases when support is insufficient, rather than universally recording them as unchecked or directly declaring normal/malignant status.

Distinguish full native depth from row sums after the declared gene projection. When discrepancies are found, first verify
omitted original-column sparse values and selection policies; neither immediately claim UID misalignment nor relax checks to conceal
unknown discrepancies. Report input columns, genomic common genes, expression-selected
genes, final feature positions, unique symbols, and repeated-position policies separately. Do not silently deduplicate
and continue using the old method name or measurement scale.

Direct differences between continuous windows require matching chromosome arms and ordered gene lists; identical
window names are insufficient for alignment. Separately report raw-arm and thresholded/trimmed-arm comparisons,
gene coverage per arm, and low-coverage limitations; unmeasured arms remain NA. When threshold summaries are reference-sensitive, do not
use them to flip normal/malignant identities. Reference eligibility for continuous RNA-CNV sensitivity analysis does not automatically extend to spatial references,
final normal identities, or DNA evidence. Record time/resources for input preparation, algorithms, export, and verification separately.

### 2.2 When normal and tumor cells can be separated

```text
tumor-like programme + concordant CNV/context evidence -> reviewed malignant support
coherent normal lineage + normal-atlas and tissue-context support -> normal-supported
CNV neutral alone -> no automatic normal label
```

### 2.3 When normal and tumor cells cannot be separated

```text
same broad epithelial/parenchymal cluster
  -> CNV abnormal arm-level pattern: supporting evidence, adjudicate with identity/reference/context
  -> coherent normal lineage plus normal-reference and tissue-context support: normal parenchymal (CNV neutrality alone is insufficient)
  -> CNV low/diploid + concordant tumor programme, patient-specific structure and tumor context: malignant-probable remains eligible
  -> genuinely conflicting or insufficient multi-source evidence: ambiguous
```

---

## 3. malignant evidence hierarchy

| evidence level | evidence type | interpretation |
|---|---|---|
| B | expression-inferred CNV / CNA from FastCNVpy, inferCNV, CopyKAT, Numbat or equivalent | supporting RNA evidence; sensitive to reference, depth and state, not independent DNA confirmation |
| A | matched WES/WGS/panel/scDNA/FISH/IHC mutation or amplification evidence | strongest external evidence |
| B | patient-specific tumor cluster plus poor normal-reference mapping | supportive |
| B | cancer-specific malignant program plus tumor-region enrichment | supportive |
| C | only canonical epithelial/parenchymal markers | not enough |
| C | only proliferation or stress | not enough |

---

## 4. malignancy_status labels

```yaml
malignant-high:
  required: convergent tumor evidence and high-confidence review; RNA-CNV alone is insufficient and high confidence is not DNA confirmation
malignant-probable:
  required: convergent tumor-context, patient-specific structure and tumor-program/normal-reference evidence; weak or absent CNV does not exclude this label
normal:
  required: coherent normal lineage plus matched-normal-reference and tissue/context support; CNV neutrality alone is insufficient
ambiguous:
  triggered_by: adequate reviewed joint evidence remains conflicting or mixed; weak CNV alone is insufficient
not_applicable:
  used_for: immune, stromal, endothelial, mural cells unless tumor lineage is expected
```

Missing evidence is `evidence_pending` workflow progress: actively retrieve or
generate it, document inaccessible evidence and do not submit a final identity or
malignancy decision. This is not a new runtime annotation-schema enum.

---

## 5. epithelial / parenchymal broad gate

| tissue context | candidate cell markers | cancer skill |
|---|---|---|
| breast | `EPCAM, TACSTD2, KRT8, KRT18, KRT19, CDH1, MUC1` | `tumor-BRCA.md` |
| colorectal | `EPCAM, TACSTD2, KRT8, KRT18, KRT19, KRT20, CDH1, CLDN3, CLDN4` | `tumor-CRC.md` |
| GBM / brain | `GFAP/AQP4/ALDOC` astrocytic, `PDGFRA/OLIG1/OLIG2` OPC-like, `SOX2/SOX4/ASCL1/DCX` NPC-like, `CD44/CHI3L1/VIM` MES-like | `tumor-GBM.md` |

---

## 6. do_not_call

```text
Do not call malignant from EPCAM/KRT/GFAP/PDGFRA alone.
Do not call malignant from RNA-CNV alone, and do not let CNV override a coherent non-tumour lineage programme.
Do not call tumour-derived endothelium, pericytes or fibroblasts (vasculogenic mimicry, tumour-to-stromal transdifferentiation) from CNV, co-clustering or shared ECM genes.
Do not apply a whole-cluster CNV call to a cluster that mixes a non-tumour lineage with tumour cells; split per cell on lineage markers.
Do not call normal from CNV-negative alone if tumor type can be CNV-low or MSI-like; use SNV/reference/tumor program.
Keep differentiation and functional states available while malignant_status is assessed separately; unresolved malignancy must not erase a supported epithelial identity.
Do not turn RNA CNV thresholds into final malignant/normal labels. Retain continuous CNV/arm evidence separately; inspect reference, low-RNA, cycling and stress effects, and resolve malignancy using the complete biological evidence.
Do not let cycling/ribosomal/stress clusters become tumor subtypes.
```

---

## 7. recommended output

```yaml
tumor_status_record:
  cluster_id: <cluster>
  broad_type: epithelial_or_parenchymal
  cnv_status: abnormal | neutral | weak | not_run
  snv_status: positive | negative | unknown
  normal_reference_mapping: matched | poor | unknown
  malignant_status: malignant-high | malignant-probable | normal | ambiguous
  next_skill: tumor-BRCA.md | tumor-CRC.md | tumor-GBM.md | normal tissue annotation
  confidence: high | medium | low
```


### Genomic axes, reference execution and algorithm semantics

CNV needs an unambiguous positional axis. Resolve duplicated symbols with verified native stable gene IDs and a build-compatible canonical position table; neither taking the first symbol column nor summing different unknown loci is defensible. Preserve unlocated native feature IDs/counts and the crosswalk. Keep original annotation-version uncertainty explicit when matching versionless stable IDs.

Save requested and actually effective reference UIDs/role labels per sample. A wrapper status contract must accept the status emitted by the real annotation committer; do not rewrite source statuses to manufacture eligibility. The immutable managed release may lag a tested canonical fix, in which case use a recorded validated project adapter. A pooled T/myeloid fit can weight role-specific factors differently from a single merged role, so do not silently change the roles.

Inspect the executed engine: FastCNV can center log2 counts on its expression-selected genomic projection rather than full-library CP10k. Export native depth, genomic counts, unique-selected counts and repeated selected-position counts separately. Ordered selected positions may contain intentional upstream-compatible repeats; they are not the same number of independent genes. Window names alone do not establish comparable ordered genes between legacy and current fits. RNA-CNV evidence and reference sensitivity remain distinct from normal/malignant assignment.
