# Overview of the Single-Cell Annotation Markers Skill

Version: v1.0-refined  
Scope: single-cell transcriptomic annotation in the current project; this does not cover complete experimental standards for spatial, protein, or multi-omics validation.  
Core purpose: files under `markers/` are neither automatic classifiers nor gene-scoring scripts; they provide an “evidence dictionary + decision workflow” for annotation subagents.

---

## 1. Goals of this refactoring

This directory retains only three categories of information:

1. **Historically established, publication-appropriate cell phenotypes or nomenclatures**  
   For example, `T_CD8_TEMRA`、`T_CD8_TEXp`、`Mac_TAM_SPP1_TREM2`、`CRC_MalignantEpi_iCMS2`、`GBM_Malignant_AC-like`、`BRCA_MalignantEpi_LuminalA`。

2. **Marker / signature evidence for these phenotypes observable in single-cell transcriptomes**  
   Use markers as a “joint evidence set”; no single gene should determine a name.

3. **The decision sequence from broad to fine annotation**  
   Each skill file includes `decision_process`, `required_evidence`, `do_not_call`, and final naming recommendations.

Functional interpretation and continua such as cell cycle, ribosomal activity, hypoxia, IFN, stress, and metabolism are in principle provided by cNMF / FastCNMF. The markers skill only maps cNMF programmes to **composable state suffixes**; it no longer designs scoring algorithms within the skill.

---

## 2. File structure

| File | Purpose |
|---|---|
| `annotation_workflow.md` | Project-level annotation workflow: DEG, cNMF, FastCNVpy, main-process adjudication and subprocess evidence preparation, and final exports |
| `naming_strategy.md` | General naming strategy: genus + differentia + state suffix + evidence level |
| `major_cell_types.md` | Broad-group annotation markers and exclusions |
| `T.md` | Historical T-cell phenotypes, markers, and 2025 modular nomenclature decision rules |
| `myeloid.md` | Monocyte/macrophage/DC/granulocyte/mast-cell annotation rules |
| `fibroblast.md` | Fibroblast / CAF annotation rules |
| `tumor-general.md` | General CNV and malignant/normal adjudication rules for tumor samples |
| `tumor-BRCA.md` | Annotation of malignant breast-cancer epithelium and normal breast parenchymal epithelium |
| `tumor-CRC.md` | Annotation of malignant colorectal-cancer epithelium and normal intestinal epithelium |
| `tumor-GBM.md` | Annotation of malignant GBM glioma cells and normal brain parenchymal cells |
| `FINAL_NAMING_REFERENCE.md` | Final naming reference table for direct consultation |

---

## 3. Global principles

### 3.1 Identity first, function second

- `cell identity`：assessed jointly using markers, exclusions, reference mapping, CNV/SNV, and tissue context.
- `functional state`：interpreted using cNMF / FastCNMF, DEG, functional enrichment, and context; it must not overwrite identity.
- `technical state`：cycling, stress, ribosomal, mitochondrial, ambient RNA, and doublet risks serve only as supplementary notes.

### 3.2 Broad-group and subgroup annotation differ

- Broad-group annotation may rely primarily on markers and exclusions.
- Fine subgroup annotation must jointly use:
  - Subgroup DEG;
  - cNMF / FastCNMF programmes;
  - Historical phenotypes in the reference literature;
  - Sample, tissue, pathology, and CNV/SNV evidence;
  - When necessary, merge similar small groups into the same canonical broad group or historical phenotype.

### 3.3 Resolve malignant / normal status first in tumor tissue

In tumor samples, parenchymal and tumor cells may be inseparable. Proceed in this order:

1. Complete broad-group annotation first.
2. Run FastCNVpy on parenchymal/epithelial/glial or other cells suspected of tumor origin.
3. If normal parenchymal and tumor cells are already separated: use FastCNVpy for supporting confirmation.
4. If they are not separated: assess malignancy jointly using continuous FastCNVpy evidence, normal atlases, tumor programmes, and patient/tissue context; low/neutral CNV must not automatically assign normal status.
5. Only malignant cells enter the cancer-specific tumor skill.
6. Name normal parenchymal cells according to normal cell lineages in the corresponding tissue.

### 3.4 Report uncertainty as uncertainty

When evidence is missing, record `evidence_pending` work progress and actively retrieve sources, tissue markers, and literature or run necessary supplementary analyses; record why evidence is unmeasurable or temporarily inexecutable, continue feasible work, and do not submit a final type. The uncertain names below are only for competing or mixed identities that remain after sufficient evidence has been collected and reviewed:

```text
Unknown_<broad_type>
Ambiguous_<broad_type>_<reason>
<cell_type>_like_lowConfidence
<cell_type>_<historical_phenotype>_mixed
```

Do not force a fine cell-subtype name onto a cluster with “only one high marker.”

---

## 4. Recommended output schema

```yaml
annotation_record:
  cluster_id: <cluster>
  broad_type: <major lineage>
  parent_type: <canonical parent cell type>
  historical_phenotype: <published phenotype if supported>
  functional_overlay: <cNMF-derived state, optional>
  malignancy_status: malignant | normal | probable_malignant | ambiguous | not_applicable
  final_label: <canonical final name>
  display_label: <short label for figures>
  positive_evidence: [marker genes, protein markers, CNV/SNV, reference labels]
  exclusion_evidence: [markers excluding adjacent lineages]
  cnmf_evidence: [program ids and top genes, if used]
  confidence: high | medium | low
  uncertainty: <remaining alternatives or reason>
  do_not_call: [labels explicitly rejected]
```
