# Fibroblast / CAF 注释 skill

版本：v1.0-refined  
范围：human tumor and normal tissue scRNA-seq。  
定位：提供 fibroblast / CAF 的真实表型与 marker；功能模块由 cNMF / FastCNMF 附加。

---

## 1. 决策流程

```text
1. Fibroblast gate: COL1A1/COL1A2/COL3A1/DCN/LUM/COL6A1/DPT/PDGFRA
2. Exclude immune, epithelial, endothelial, pericyte/SMC
3. Decide resident fibroblast vs activated CAF-like
4. Assign historical phenotype: resident PI16/DPT, matrix CAF, myCAF, iCAF, apCAF, LRRC15+ CAF, reticular/TLS CAF, etc.
5. Attach cNMF function/niche overlays: ECM, contractile, inflammatory, IFN, hypoxia, cycling, TLS, vascular
6. For cancer-specific interpretation, add CRC/BRCA/PDAC context when evidence supports it
```

---

## 2. Fibroblast gate and exclusions

| decision | markers |
|---|---|
| fibroblast lineage positive | `COL1A1, COL1A2, COL3A1, DCN, LUM, COL6A1, DPT, PDGFRA, PDPN` |
| normal/resident-like | `PI16, DPP4, CD34, MFAP5, COL15A1, DPT` |
| activated CAF-like | `FAP, PDPN, POSTN, MMP11, COL11A1, CTHRC1, LRRC15` |
| exclude immune | `PTPRC, LST1, CD68, CD3D, MS4A1` |
| exclude epithelial/tumor | `EPCAM, KRT8, KRT18, KRT19, CDH1` |
| exclude endothelial | `PECAM1, VWF, KDR, CLDN5` |
| exclude pericyte/SMC | `RGS5, CSPG4, PDGFRB, MYH11, MCAM, ACTA2/TAGLN` if mural program dominates |

`FAP`、`ACTA2`、`PDPN` 不可单独定义 CAF subtype。

---

## 3. Resident fibroblast phenotypes

| phenotype | markers | final label | notes |
|---|---|---|---|
| PI16/DPP4 resident fibroblast | `PI16, DPP4, CD34, MFAP5, CFD, C3, COL15A1` | `FB_Resident_PI16_DPP4` | common interstitial/resident-like fibroblast |
| DPT/LUM matrix resident fibroblast | `DPT, LUM, DCN, COL1A1, COL1A2, COL6A1, PDGFRA` | `FB_Resident_DPT_LUM` | basal fibroblast identity |
| COL15A1 vascular-associated fibroblast | `COL15A1, PI16, MFAP5, CD34`, low endothelial markers | `FB_Resident_COL15A1_VascularAssoc` | distinguish from endothelial/pericyte |
| Cycling fibroblast | parent fibroblast markers + `MKI67, TOP2A, UBE2C` | `FB_<parent>_Cycling` | state only |

---

## 4. CAF phenotypes

| phenotype | markers | final label | notes |
|---|---|---|---|
| Matrix / ECM-remodeling CAF | `MMP11, CDH11, POSTN, COL1A1, COL1A2, COL1A3, COL11A1, FN1, TNC` | `CAF_Matrix_POSTN_COL11A1` | ECM deposition/remodeling, invasive stroma |
| Contractile myCAF | `ACTA2, TAGLN, MYL9, TPM1, TPM2, CNN1, CTGF` | `CAF_myCAF_Contractile` | if `RGS5/CSPG4/PDGFRB/MYH11` dominates, call mural instead |
| LRRC15+ TGFβ-driven CAF | `LRRC15, POSTN, THBS2, COL11A1, CTHRC1` | `CAF_LRRC15_TGFb_Matrix` | immunosuppressive matrix-associated CAF phenotype |
| Inflammatory iCAF | `IL6, LIF, CXCL1, CXCL2, CXCL8, CXCL12, CCL2, PLA2G2A, CFD, C3, CD34` | `CAF_iCAF_CXCL12_IL6` | cytokine/chemokine CAF; distinguish ambient inflammation |
| Complement/secretory CAF | `C3, C7, CFB, CFD, CFH, CFI` | `CAF_Complement_Inflammatory` | can overlap with iCAF |
| Antigen-presenting CAF | `CD74, HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, CIITA`; must be `PTPRC-, LST1-, CD68-` | `CAF_apCAF_MHCII` | exclude myeloid/B-cell doublet |
| IFN-response CAF | `IDO1, ISG15, MX1, OAS1, IFIT1, IFIT2, IFIT3, CXCL9, CXCL10, CXCL11` | `CAF_IFNhi` | state overlay |
| Vascular-support CAF | `MCAM, COL18A1, NOTCH3, NID2`, low/int `ACTA2` | `CAF_VascularSupport` | if `RGS5/CSPG4/PDGFRB` high, call pericyte |
| Reticular / TLS CAF | `CCL19, CCL21, CXCL13, IL7, PDPN, LTBR` | `CAF_TLS_Reticular` | stronger with TLS/B-cell follicle context |
| Tumor-border / hypoxic CAF | `MME/CD10, NT5E/CD73, ENO1, NDRG1, CA9, HSPH1, PDPN` | `CAF_TumorBorder_Hypoxic` | strongest with spatial tumor-border evidence |
| Cycling/stress CAF | `MKI67, TOP2A, UBE2C`; `HSPA1A, HSPH1, DNAJB1` | `CAF_<parent>_Cycling` or `CAF_<parent>_Stress` | state only |

---

## 5. Cancer-context rules

### 5.1 CRC

CRC fibroblast/CAF labels should preserve crypt niche and WNT/BMP axes when present:

| phenotype | markers | final label |
|---|---|---|
| CRC matrix CAF / CAF-A-like | `MMP2, DCN, COL1A2, COL11A1, POSTN` | `CRC_CAF_Matrix_ECM` |
| CRC contractile CAF / CAF-B-like | `ACTA2, TAGLN, PDGFA, MYL9` | `CRC_CAF_myCAF_Contractile` |
| Crypt-bottom fibroblast | `WNT2, WNT2B, RSPO3, GREM1`, often `PDGFRAlow` | `CRC_FB_CryptBottom_WNT` |
| Crypt-top fibroblast | `WNT5A`, BMP ligands, often `PDGFRAhigh` | `CRC_FB_CryptTop_BMP` |
| Tumor-border hypoxic CAF | `CA9, NDRG1, ENO1, MME, NT5E` | `CRC_CAF_TumorBorder_Hypoxic` |

### 5.2 BRCA

| phenotype | markers | final label |
|---|---|---|
| Perivascular CAF | `NID2, MCAM, NOTCH3, COL18A1` | `BRCA_CAF_Perivascular` |
| Collagen-rich matrix CAF | `PDGFRA, MMP11, POSTN, COL1A1, COL11A1` | `BRCA_CAF_Matrix` |
| Immune-suppressive CAF-S1-like | `FAP, PDPN, CXCL12, IL6, PDGFRB, DPP4, NT5E` | `BRCA_CAF_ImmuneSuppressive_CXCL12` |
| TLS-adjacent reticular CAF | `CCL19, CCL21, CXCL13, PDPN` | `BRCA_CAF_TLS_Reticular` |
| Tumor-border hypoxic CAF | `MME, NT5E, NDRG1, ENO1, CA9` | `BRCA_CAF_TumorBorder_Hypoxic` |

### 5.3 PDAC reference pattern

PDAC often uses:

```text
PDAC_CAF_myCAF_Contractile        # ACTA2/TAGLN/CTGF, tumor-proximal
PDAC_CAF_iCAF_CXCL12_IL6          # IL6/LIF/CXCLs, inflammatory
PDAC_CAF_apCAF_MHCII              # CD74/MHC-II, PTPRC-negative
PDAC_CAF_LRRC15_TGFb_Matrix       # LRRC15/POSTN/THBS2/COL11A1/CTHRC1
```

---

## 6. Conflict rules

```text
CD74/HLA-DRA high + PTPRC/LST1/CD68 high -> myeloid/B doublet or immune cell, not apCAF.
ACTA2/TAGLN high + RGS5/CSPG4/PDGFRB/MYH11 high -> pericyte/SMC before myCAF.
COL1A1/COL1A2 high + EPCAM/KRT high + CNV abnormal -> EMT tumor possible, not CAF by default.
FAP high alone -> activated fibroblast/CAF-like only, not subtype.
MKI67/TOP2A high -> cycling overlay only.
HSPA1A/FOS/JUN high across all lineages -> dissociation stress.
```

---

## 7. Final labels

```text
FB_Resident_PI16_DPP4
FB_Resident_DPT_LUM
FB_Resident_COL15A1_VascularAssoc
CAF_Matrix_POSTN_COL11A1
CAF_myCAF_Contractile
CAF_LRRC15_TGFb_Matrix
CAF_iCAF_CXCL12_IL6
CAF_Complement_Inflammatory
CAF_apCAF_MHCII
CAF_IFNhi
CAF_TLS_Reticular
CAF_TumorBorder_Hypoxic
CAF_Unresolved
```
