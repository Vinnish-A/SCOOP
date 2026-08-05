# 髓系细胞注释 skill：单核、巨噬、DC、粒细胞、肥大细胞

版本：v1.0-refined  
范围：human tumor / normal tissue scRNA-seq；可作为 mouse 数据的 marker 起点但需物种复核。  
定位：提供真实历史表型与 marker；功能连续谱由 cNMF / FastCNMF 附加。

---

## 1. 决策流程

```text
1. Myeloid gate: LYZ/LST1/TYROBP/FCER1G/AIF1/CST3/CTSS
2. Exclude T/B/NK/epithelial/stromal/endothelial doublets
3. Split into monocyte / macrophage / DC / neutrophil / mast / ambiguous myeloid
4. For macrophage: separate origin-like evidence from state evidence
5. For tumor macrophage: report TAM state but do not infer absolute ontogeny without lineage evidence
6. Add cNMF overlays: IFN, APC, SPP1/TREM2, cycling, stress, hypoxia, lipid/phagocytic
7. Never use M1/M2 as final labels
```

---

## 2. Myeloid gate

| category | markers |
|---|---|
| pan-myeloid | `LYZ, LST1, TYROBP, FCER1G, AIF1, CST3, CTSS, SPI1` |
| macrophage core | `C1QA, C1QB, C1QC, APOE, MERTK, MSR1, CD68, CD163` |
| monocyte core | `FCN1, VCAN, S100A8, S100A9, CCR2, LILRB1` |
| DC core | `FCER1A, CD1C, CLEC10A, CLEC9A, XCR1, LAMP3, CCR7, IRF8` |
| exclude T | `CD3D, CD3E, TRAC` |
| exclude B/plasma | `MS4A1, CD79A, CD79B, MZB1, XBP1` |
| exclude epithelial | `EPCAM, KRT8, KRT18, KRT19` |
| exclude fibroblast | `COL1A1, DCN, LUM, PDGFRA` |

---

## 3. Monocyte phenotypes

| phenotype | markers | final label | notes |
|---|---|---|---|
| Classical / inflammatory monocyte | `FCN1, VCAN, S100A8, S100A9, S100A12, CCR2, IL1B, CXCL8, LILRB1` | `Mono_FCN1_Classical` or `Mac_Mo_FCN1_IL1B` if macrophage-like | common entry/inflammatory axis in tumors |
| Non-classical monocyte | `FCGR3A, MS4A7, CX3CR1, LST1, LILRB2, AIF1` | `Mono_FCGR3A_Nonclassical` | check NK contamination if `NKG7/GNLY` high |
| Intermediate monocyte | `FCN1` with `FCGR3A/MS4A7`, MHC-II genes | `Mono_Intermediate` | often transitional; do not over-interpret |
| Cycling monocyte | parent monocyte markers + `MKI67, TOP2A, UBE2C` | `Mono_<parent>_Cycling` | state only |

---

## 4. Macrophage origin-like and tissue programs

| phenotype | markers | final label | notes |
|---|---|---|---|
| C1QC phagocytic macrophage | `C1QA, C1QB, C1QC, APOE, MERTK, MSR1, CD68, CTSB, CTSD` | `Mac_C1QC_Phagocytic` | broad macrophage state; can be resident-like or TAM-like |
| Mo-derived-like macrophage | `FCN1, VCAN, S100A8, S100A9, CCR2, IL1B, CXCL8, LILRB1` plus macrophage core | `Mac_Mo_FCN1_IL1B` | origin-like, not proof of ontogeny |
| TRM-like / resident-like macrophage | `FOLR2, LYVE1, MRC1, TIMD4, C1QC, SEPP1, GAS6, VSIG4, F13A1, MARCO` | `Mac_TRM_FOLR2_LYVE1` | tissue-resident-like; marker strength varies by organ |
| Microglia-like TRM | `P2RY12, TMEM119, SALL1, GPR34, CX3CR1, OLFML3, CSF1R, AIF1` | `Mac_Microglia_P2RY12` | GBM/brain; activated microglia may lose homeostatic markers |
| SPP1/TREM2 TAM | `SPP1, TREM2, APOE, GPNMB, LGALS3, FABP5, CTSB, CTSD` | `Mac_TAM_SPP1_TREM2` | tumor-educated remodeling/lipid/phagocytic state; not absolute origin |
| APC-like macrophage | `HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, CD74, CIITA` | `Mac_APC_MHCII` | distinguish from cDC and B cells |
| IFN-response macrophage | `IFIT1, IFIT2, IFIT3, ISG15, IFI6, IFI44L, MX1, OAS1` | `Mac_IFNhi` | overlay/state only |
| Cycling macrophage | parent macrophage markers + `MKI67, TOP2A, STMN1, TYMS` | `Mac_<parent>_Cycling` | overlay/state only |

---

## 5. DC phenotypes

| phenotype | markers | final label | do_not_call |
|---|---|---|---|
| cDC1 | `CLEC9A, XCR1, BATF3, IRF8, CADM1, THBD` | `DC_cDC1` | do not call from `IRF8` alone |
| cDC2 | `CD1C, FCER1A, CLEC10A, IRF4, CD1E, HLA-DRA` | `DC_cDC2` | distinguish from monocytes with `FCN1/VCAN` |
| Migratory / mature DC | `LAMP3, CCR7, CCL19, FSCN1, CD83, CD40, RELB` | `DC_LAMP3_Migratory` | state may derive from cDC1/cDC2; use as phenotype |
| pDC | `LILRA4, IL3RA, CLEC4C, TCF4, IRF7, GZMB` | `DC_pDC` | high `GZMB` alone not enough |
| moDC-like | `FCN1, VCAN, CD1C, FCER1A, CLEC10A, LYZ` | `DC_moDC_like` | low confidence unless reference supports |

---

## 6. Neutrophil / mast / granulocyte-like cells

| phenotype | markers | final label | notes |
|---|---|---|---|
| Neutrophil | `S100A8, S100A9, CSF3R, CXCR2, FCGR3B, MPO, ELANE, CEACAM8` | `Neutrophil` | may be low-quality / dissociation sensitive |
| Low-density granulocyte / PMN-MDSC-like | `S100A8, S100A9, ARG1, FCGR3B, CXCR2, LCN2` | `Neutrophil_PMNMDSC_like` | only use in tumor/inflammation context; low confidence without function |
| Mast cell | `KIT, TPSAB1, TPSB2, CPA3, MS4A2, HDC` | `Mast` | exclude basophil if marker set available |
| Basophil-like | `FCER1A, MS4A2, CPA3, HDC, IL3RA` without strong `KIT/TPSAB1` | `Basophil_like` | rare in tissue scRNA; review |

---

## 7. Tumor macrophage minimum labels by cancer context

| cancer context | minimum macrophage labels | marker focus |
|---|---|---|
| GBM | `Mac_Microglia_P2RY12`, `Mac_Mo_FCN1_IL1B`, `Mac_TAM_SPP1_TREM2` | `P2RY12/TMEM119/SALL1`; `FCN1/S100A8/S100A9/CCR2`; `SPP1/APOE/TREM2/GPNMB` |
| BRCA | `Mac_Mo_FCN1_IL1B`, `Mac_C1QC_FOLR2`, `Mac_CD274_PDCD1LG2`, `Mac_TAM_SPP1_TREM2` | `FCN1/S100A8`; `C1QC/FOLR2`; `CD274/PDCD1LG2`; `SPP1/TREM2` |
| CRC | `Mac_Inflammatory_IL1B`, `Mac_Proangiogenic_Remodeling`, `Mac_APC_MHCII` | `IL1B/CXCL8`; `SPP1/MMP9/VEGFA`; `HLA-DRA/CD74` |
| HCC | `Mac_C1QC_Phagocytic`, `Mac_MMP9_TerminalTAM`, `Mac_TAM_SPP1_TREM2` | `C1QC`; `MMP9`; `SPP1/TREM2` |
| PDAC | `Mac_Mo_TumorContinuum`, `Mac_TAM_SPP1_TREM2`, `Mac_C1QC_ResidentLike` | compare tumor and PBMC when available |

---

## 8. Final naming examples

```text
Mono_FCN1_Classical
Mono_FCGR3A_Nonclassical
Mac_Mo_FCN1_IL1B
Mac_C1QC_Phagocytic
Mac_TRM_FOLR2_LYVE1
Mac_Microglia_P2RY12
Mac_TAM_SPP1_TREM2
Mac_APC_MHCII
Mac_IFNhi
DC_cDC1
DC_cDC2
DC_LAMP3_Migratory
DC_pDC
Mast
Neutrophil
```

---

## 9. do_not_call

```text
Do not use M1/M2 as final macrophage labels.
Do not infer monocyte-derived vs tissue-resident absolute origin without lineage, fate, genotype or strong reference evidence.
Do not call macrophage from LYZ alone.
Do not call DC from CD74/HLA-DRA alone.
Do not call microglia absent because P2RY12/TMEM119 are down; activated microglia may lose homeostatic markers.
Do not call SPP1/TREM2 as an origin; it is a tumor-education/remodeling state.
Do not create a new cell type for IFN, cycling, ribosomal or stress programs.
```
