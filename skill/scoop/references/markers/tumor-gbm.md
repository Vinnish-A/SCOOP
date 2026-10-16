# GBM 肿瘤细胞注释 skill：恶性胶质瘤细胞与正常脑实质细胞

版本：v1.0-refined  
范围：human IDH-wildtype GBM / glioma-like scRNA/snRNA-seq。  
定位：必须先区分正常脑实质细胞与 malignant glioma cells，再做 AC/OPC/NPC/MES 等状态注释。

---

## 1. 决策流程

```text
1. Broad brain/parenchymal annotation: astrocyte, OPC, oligodendrocyte, neuron, ependymal/ciliated, microglia, endothelial, mural, fibroblast/VLMC
2. Candidate glioma cells: marker overlap with AC-like, OPC-like, NPC-like, MES-like programs
3. FastCNVpy / SNV / tumor reference mapping to assign malignant_status
4. CNV-low/neutral cells -> joint review of genomic, normal-atlas, tumour-programme and patient/tissue evidence; no automatic normal routing
5. malignant cells -> GBM historical states: AC-like, OPC-like, NPC-like, MES-like
6. Reference-mapping may support GBmap-like states, but final identity requires model/analyst review rather than classifier naming
7. cNMF programs -> overlays: cycling, hypoxia, stress, gliosis, ECM, connectivity, invasion, cilia-like, neuronal-like
```

---

## 2. Malignancy evidence

GBM marker overlap is severe; malignant status must prioritize genomic evidence.

| evidence | markers / events | interpretation |
|---|---|---|
| CNV/CNA | `chr7 gain`, `chr10 loss`, `EGFR amplification`, `CDKN2A/B deletion`, `PDGFRA amplification`, `CDK4/MDM2 amplification`, `PTEN region loss` | strongest scRNA-level evidence |
| SNV/fusion | tumor WES/WGS/panel, EGFRvIII, allele-specific CNA | strong evidence |
| tumor reference mapping | GBM malignant-state reference such as AC/OPC/NPC/MES | supportive; resolve normal-cell overlap with joint genomic, patient/tissue and programme evidence, not a mandatory positive CNV gate |
| marker alone | `GFAP`, `PDGFRA`, `SOX2`, `CD44`, etc. | insufficient |

---

## 3. Normal brain parenchymal and TME cells

Normal origin requires coherent brain lineage plus normal-reference and tissue-context support. CNV neutrality or absent/unassayed tumour SNVs alone does not establish normality.

| normal cell | markers | final label | confusion |
|---|---|---|---|
| Astrocyte | `AQP4, ALDH1L1, SLC1A2, SLC1A3, GLUL, GJA1, GFAP, S100B, MLC1` | `GBM_Normal_Astrocyte` | overlaps AC-like tumor |
| OPC | `PDGFRA, CSPG4, VCAN, PTPRZ1, OLIG1, OLIG2, SOX10` | `GBM_Normal_OPC` | overlaps OPC-like tumor |
| Oligodendrocyte | `MBP, PLP1, MOG, MAG, MOBP, OPALIN, CLDN11` | `GBM_Normal_Oligodendrocyte` | myelin RNA contamination possible |
| Neuron | pan-neuronal `RBFOX3, SNAP25, SYT1, MAP2, TUBB3`; excitatory `SLC17A7, CAMK2A`; inhibitory `GAD1, GAD2, SLC6A1` | `GBM_Normal_Neuron` | overlaps NPC/NEU-like tumor/ambient RNA |
| Ependymal / ciliated | `FOXJ1, PIFO, TPPP3, DNAH5, CFAP*` | `GBM_Normal_EpendymalCiliated` | overlaps malignant cilia-like state if CNA+ |
| Microglia | `P2RY12, TMEM119, CX3CR1, SALL1, CSF1R, AIF1` | `Mac_Microglia_P2RY12` | myeloid skill, not tumor |
| Endothelial | `PECAM1, VWF, CLDN5, KDR, FLT1` | `Endo_Brain` | vascular niche |
| Pericyte / mural | `PDGFRB, RGS5, MCAM, ACTA2, TAGLN` | `Pericyte_Brain` or `SMC_Brain` | overlaps MES/ECM tumor if CNA+ |
| Fibroblast / VLMC | `COL1A1, COL1A2, DCN, LUM, COL3A1` | `FB_VLMC` | overlaps MES-ECM tumor if CNA+ |

---

## 4. GBM malignant historical states

Only use these labels for `malignant-high` or `malignant-probable` cells.

### 4.1 AC-like / astrocyte-like

| subtype / overlay | markers | final label | notes |
|---|---|---|---|
| AC-like main state | `GFAP, AQP4, ALDOC, APOE, SLC1A3, S100B, MLC1, SOX9` | `GBM_Malignant_AC-like` | CNA/SNV required because normal astrocytes overlap |
| AC-homeostatic / differentiated | `AQP4, ALDOC, SLC1A2, SLC1A3, GLUL, GFAP`; low `MKI67/TOP2A` | `GBM_Malignant_AC-like_Homeostatic` | differentiated-like overlay |
| AC-gliosis / injury-response | `VIM, SERPINA3, LCN2, C3, CD44, CHI3L1, GFAP` | `GBM_Malignant_AC-like_Gliosis` | injury/gliosis overlay |
| AC-connectivity / tumor-microtube-like | `CHI3L1, GAP43, APOE, GJA1, TTYH1` | `GBM_Malignant_AC-like_Connectivity` | stronger with spatial/morphology evidence |
| Cilia-like | `DNAH*, CFAP*, FOXJ1, TPPP3` | `GBM_Malignant_AC-like_CiliaLike` | require joint tumor evidence; CNA-negative alone does not establish normal ependymal origin |

### 4.2 OPC-like

| subtype / overlay | markers | final label | notes |
|---|---|---|---|
| OPC-like main state | `PDGFRA, OLIG1, OLIG2, SOX10, CSPG4, PTPRZ1, BCAN, VCAN, NKX2-2` | `GBM_Malignant_OPC-like` | PDGFRA amplification supports |
| OPC-progenitor / stem-like | `PDGFRA, CSPG4, PTPRZ1, OLIG2, SOX10, NES, PROM1` | `GBM_Malignant_OPC-like_Progenitor` | progenitor overlay |
| OPC-cycling | `MKI67, TOP2A, PCNA, MCM2-7, CENPF, UBE2C, CDK1, CCNB1, CCNB2` | `GBM_Malignant_OPC-like_Cycling` | cycling is overlay, not fifth classic state |
| OPC-invasive / motile | `TNC, ANXA1, ITGA6, ITGB1, MMP2, MMP14, SERPINE1` | `GBM_Malignant_OPC-like_Invasive` | invasion overlay |
| Oligo-differentiating tumor | `PLP1, MBP, MAG, MOG, MOBP, OPALIN` | `GBM_Malignant_OPC-like_OligoDifferentiating` | require joint tumor support beyond differentiation markers; CNA-low remains eligible |

### 4.3 NPC-like / neural progenitor-like

| subtype / overlay | markers | final label | notes |
|---|---|---|---|
| NPC-like main state | `SOX2, SOX4, SOX11, ASCL1, DLL3, HES6, DCX, STMN2, TUBB3, NEUROD1` | `GBM_Malignant_NPC-like` | CDK4 amplification supports |
| NPC1-like developmental/progenitor | `SOX2, ASCL1, DLL3, HES6, SOX4, SOX11` | `GBM_Malignant_NPC1-like` | use if reference mapping supports |
| NPC2-like neuronal-lineage | `DCX, STMN2, TUBB3, MAP2, NEUROD1` | `GBM_Malignant_NPC2-like` | distinguish normal neuron |
| NEU-like updated state | `NRG1, NRG3, NRXN3, CNTNAP2, SYN3, SYT1` | `GBM_Malignant_NEU-like` | use as secondary/reference-mapped label unless strong classifier support |
| NPC-cycling | `MKI67, TOP2A, MCM2-7, UBE2C, BIRC5` | `GBM_Malignant_NPC-like_Cycling` | overlay |

### 4.4 MES-like / mesenchymal-like

| subtype / overlay | markers | final label | notes |
|---|---|---|---|
| MES-like main state | `CD44, CHI3L1, VIM, LGALS3, S100A10, ANXA1, SERPINE1, TNC, ITGA5, STAT3, CEBPB` | `GBM_Malignant_MES-like` | NF1 loss/mutation supports |
| MES1-inflammatory / gliosis-like | `CD44, CHI3L1, LGALS3, S100A10, ANXA1, VIM, CCL2, CXCL8, IL6, STAT3, CEBPB` | `GBM_Malignant_MES-like_InflammatoryGliosis` | macrophage-rich niche may support |
| MES2-hypoxia / angiogenic | `VEGFA, CA9, HILPDA, ADM, BNIP3, NDRG1, SLC2A1, LDHA, EGLN3` | `GBM_Malignant_MES-like_Hypoxia` | stronger near necrosis/hypoxic niche |
| MES-stress / UPR | `HSPA1A, HSPA1B, HSP90AA1, DNAJB1, ATF3, DDIT3, XBP1` | `GBM_Malignant_MES-like_Stress` | check dissociation artifact |
| MES-ECM / invasive-remodeling | `TNC, FN1, COL1A1, COL1A2, COL3A1, ITGA5, MMP2, MMP14, SERPINE1` | `GBM_Malignant_MES-like_ECMInvasive` | if no CNA, consider fibroblast/VLMC |
| MES-connectivity | `CHI3L1, GAP43, APOE` | `GBM_Malignant_MES-like_Connectivity` | requires supporting context |

---

## 5. Reference mapping layer

Reference models may output more detailed GBM states. Store them separately from the canonical four-state label.

```yaml
classic_state: AC-like | OPC-like | NPC-like | MES-like | hybrid | low_confidence
reference_state: AC-like | OPC-like | NPC1-like | NPC2-like | MES1-like | MES2-like | NEU-like | GPC-like | Cilia-like | Hypoxia | Stress | unknown
functional_overlay: cycling | hypoxia | stress | gliosis | ECM | invasion | connectivity | cilia | neuronal-like
```

Rules:

```text
ClassicState can fold NPC1/NPC2 into NPC-like and MES1/MES2 into MES-like.
ReferenceState can preserve higher-resolution labels when the classifier or project reference supports them.
Cycling, hypoxia and stress remain overlays unless the reference model explicitly treats them as malignant-state labels.
Hybrid labels are allowed when top state evidence is close and doublet risk is low.
```

---

## 6. Final naming examples

```text
GBM_Normal_Astrocyte
GBM_Normal_OPC
GBM_Normal_Oligodendrocyte
GBM_Normal_Neuron
GBM_Malignant_AC-like_CNV7gain10loss
GBM_Malignant_AC-like_Gliosis
GBM_Malignant_OPC-like_Cycling
GBM_Malignant_NPC-like_Developmental
GBM_Malignant_NPC2-like_NeuronalLineage
GBM_Malignant_MES-like_InflammatoryGliosis
GBM_Malignant_MES-like_Hypoxia
GBM_Malignant_AC_MES_hybrid
GBM_PutativeMalignant_OPC-like_CNVweak
```

---

## 7. do_not_call

```text
Do not call AC-like tumor from GFAP/AQP4 alone.
Do not call OPC-like tumor from PDGFRA/OLIG2 alone.
Do not call NPC/NEU-like tumor from STMN2/SYT1/SNAP25 alone; normal neuron/ambient RNA is common.
Do not call MES-like tumor from CD44 alone.
Do not call MES-ECM tumor without CNA/SNV if COL1A1/DCN/LUM suggest fibroblast/VLMC.
Do not make cycling a primary GBM state in this project; use overlay.
Do not ignore immune-tumor doublet risk when PTPRC and tumor CNV/EGFR are both high.
```
