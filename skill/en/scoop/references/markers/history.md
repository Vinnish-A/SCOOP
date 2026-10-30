# Markers Skill Refactoring Record

Version: v1.0-refined

---

## 1. Structural changes

- Added `annotation_workflow.md`: formalizes the project-level annotation workflow, FastCNVpy branches, and the role of cNMF/FastCNMF as an SOP.
- Added `naming_strategy.md`: defines the general naming strategy of “genus + differentia + state suffix + evidence level.”
- Added `major_cell_types.md`: provides a marker dictionary and exclusions for broad-group annotation.
- Added `tumor-general.md`: abstracts malignant / normal adjudication in tumor samples into general rules.
- Added `FINAL_NAMING_REFERENCE.md`: provides a final naming table for consultation.

---

## 2. Content changes

- T cells now follow the 2025 modular nomenclature style: historical phenotype + marker + modular attributes, without over-reifying individual states.
- Myeloid removes M1/M2 as final labels, replacing them with hierarchical annotation of monocyte / macrophage / DC / TAM states.
- Fibroblast / CAF retains real CAF phenotypes and cancer-specific niches, while treating ECM, contractile, IFN, cycling, and similar signals as state or functional overlays.
- BRCA / CRC / GBM tumor cells are handled separately by cancer type; every malignant subtype first requires malignant_status.
- CNV provides malignant/normal and clonal evidence, rather than directly interpreting functional states.
- cNMF / FastCNMF is used for functional continua and identity disambiguation, not as a hard-threshold scorer within the skill.

---

## 3. Removed or downgraded content

- Removed the tendency to name directly from single-gene scores.
- Removed or downgraded statements equating “function with cell type.”
- `Cycling`, `RibosomalHigh`, `Stress`, `Hypoxiahi`, `IFNhi`, and similar labels are consistently downgraded to overlays.
- `CMS1-4` is no longer a CRC single-cell epithelial label; it is only a sample/region pseudobulk interpretation layer.
- `Normal-like` is no longer a main BRCA tumor label; by default it enters normal/ambient/doublet/CNV review.
- GBM `Cycling` is not a fifth classical state; it is an overlay on AC/OPC/NPC/MES.

---

## 4. Expected subagent behavior

```text
Read marker evidence + DEG + cNMF programs + FastCNVpy + project context
  -> Find the corresponding markers skill
  -> First assess identity and malignant_status
  -> Then map the historical phenotype
  -> Finally attach the functional overlay
  -> Output final_label, confidence, uncertainty, do_not_call
```
