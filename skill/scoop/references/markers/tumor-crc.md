# CRC / COAD / READ 肿瘤细胞注释 skill

版本：v1.0-refined  
范围：human colorectal cancer malignant epithelial and normal intestinal epithelial scRNA-seq。  
定位：只在 epithelial 候选细胞上使用；malignant_status 综合遗传、正常图谱、患者结构与组织背景判断。FastCNVpy 是支持证据，不是所有恶性细胞必须通过的阳性门槛。

---

## 1. 决策流程

```text
1. Epithelial gate: EPCAM/TACSTD2/KRT8/KRT18/KRT19/KRT20/CDH1/CLDN3/CLDN4
2. Exclude immune, fibroblast, endothelial, mural doublets
3. Assign malignant_status using FastCNVpy/SNV/reference/patient-specific clustering
4. CNV-low / neutral cells -> joint review of normal reference, tumor programme and patient/tissue context; never automatically normal
5. malignant cells -> iCMS2 / iCMS3 / mixed / indeterminate
6. Add cNMF overlays: stem-like, TA-proliferative, goblet-mucinous, IFN, metabolic, EMT-invasive, fetal-regenerative, hypoxia, DTP-like
7. CMS/CRIS/IMF/PDS are orthogonal/pseudobulk layers, not direct per-cell labels
```

---

## 2. Epithelial gate and exclusions

| item | markers |
|---|---|
| epithelial positive | `EPCAM, TACSTD2, KRT8, KRT18, KRT19, KRT20, CDH1, CLDN3, CLDN4, MSLN` |
| EMT/tumor-budding caution | `EPCAM` may decrease but `KRT8/KRT18/KRT19` remain |
| exclude immune | `PTPRC, CD3D, CD3E, TRAC, MS4A1, CD79A, LST1, FCGR3A, CD68` |
| exclude fibroblast/CAF | `COL1A1, COL1A2, DCN, LUM, COL3A1, PDGFRA, ACTA2, TAGLN` |
| exclude endothelial | `PECAM1, VWF, KDR, EMCN, PLVAP` |
| exclude mural | `RGS5, MCAM, CSPG4, MYH11` |

---

## 3. Normal intestinal epithelial phenotypes

Normal labels require coherent intestinal lineage plus matched-normal-reference and tissue-context support. CNV neutrality or an unassayed/undetected SNV is not proof of normality; CNV-low malignant candidates remain in the malignancy review.

| phenotype | markers | final label |
|---|---|---|
| Stem-like / crypt-base | `LGR5, OLFM4, ASCL2, SMOC2, AXIN2, RNF43, LRIG1` | `CRC_NormalEpi_CryptStem` |
| Transit-amplifying / proliferating | `MKI67, TOP2A, PCNA, MCM2, MCM5, TYMS, CENPF, UBE2C` | `CRC_NormalEpi_TA` |
| Absorptive colonocyte | `CA1, CA2, AQP8, SLC26A3, CEACAM7, MS4A12, GUCA2A, GUCA2B, KRT20` | `CRC_NormalEpi_AbsorptiveColonocyte` |
| Goblet | `MUC2, TFF3, SPINK4, CLCA1, FCGBP, AGR2` | `CRC_NormalEpi_Goblet` |
| BEST4+ absorptive | `BEST4, OTOP2, SPIB, GUCA2B, CA7, NPY` | `CRC_NormalEpi_BEST4Absorptive` |
| Enteroendocrine | `CHGA, CHGB, NEUROD1, PAX6, PYY, GCG, TPH1, SCG5` | `CRC_NormalEpi_Enteroendocrine` |
| Tuft | `DCLK1, TRPM5, POU2F3, AVIL, SH2D6, GFI1B` | `CRC_NormalEpi_Tuft` |
| Paneth-like / deep crypt secretory | `DEFA5, DEFA6, LYZ, REG3A, PLA2G2A, MMP7, KIT` | `CRC_NormalEpi_PanethLike_DCS` |

---

## 4. Malignant status logic

```yaml
normal_epithelial:
  epithelial_markers: high
  cnv: negative
  tumor_snv: negative_or_unknown
  reference: maps_to_normal_lineage
  patient_specific_tumor_cluster: absent
malignant_high:
  epithelial_markers: high_or_EMT_epithelial_residual
  cnv: positive_or_tumor_snv_positive
  patient_specific_tumor_cluster: present
malignant_probable:
  epithelial_markers: high
  cnv: negative_or_weak
  tumor_program: high
  normal_reference_mapping: poor
  note: iCMS3/MSI-like tumors can be CNV-low
ambiguous:
  reasons: low_UMI, doublet, weak_CNV, only_single_marker, conflicting_normal_and_tumor_programs
```

CRC CIN-supportive CNV patterns:

```text
gains: 7p/q, 8q, 13q, 20p/q
losses: 1p, 4p/q, 8p, 14q, 15q, 17p, 18p/q
```

Common driver / status evidence:

```text
APC, KRAS, NRAS, BRAF, TP53, PIK3CA, SMAD4, FBXW7, TCF7L2, ACVR2A, TGFBR2, MLH1/MSH2/MSH6/PMS2 status
```

---

## 5. Malignant epithelial intrinsic subtype

### 5.1 iCMS2

| program | markers |
|---|---|
| WNT/TCF | `ASCL2, LGR5, AXIN2, RNF43, NOTUM, NKD1, TCF7, LEF1` |
| MYC/E2F/cell cycle | `MYC, E2F1, E2F2, MKI67, TOP2A, MCM2-7, PCNA, TYMS, UBE2C` |
| EGFR ligand axis | `EREG, AREG` |
| CIN support | chr8q/13q/20q gains, chr17p/18q losses |

Final label:

```text
CRC_MalignantEpi_iCMS2
```

Common overlays:

```text
StemLike
TAProliferative
EGFRLigand
PDS1-like
PDS3-like differentiated slow-cycling
```

### 5.2 iCMS3

| program | markers |
|---|---|
| secretory/goblet | `MUC2, TFF3, SPINK4, CLCA1, FCGBP, AGR2, REG4` |
| metabolic | `SLC2A1, LDHA, HK2, PGK1, FABP1, APOA1, APOA4, MTTP, ALDOB` |
| inflammatory/IFN | `IFIT1, IFIT2, IFIT3, ISG15, MX1, OAS1, STAT1, IRF1, CXCL10, HLA-A/B/C, B2M` |
| EMT/invasive support | `LAMC2, LAMB3, ITGA5, MMP7, MMP14, SERPINE1, TGFBI, EMP1, L1CAM` |

Final label:

```text
CRC_MalignantEpi_iCMS3
```

Common overlays:

```text
GobletMucinous
InflammatoryIFN
Metabolic
EMTInvasive
FetalRegenerative
HypoxiaStress
DTPLike
```

### 5.3 mixed / indeterminate

```text
CRC_MalignantEpi_iCMSmixed
CRC_MalignantEpi_iCMSindeterminate
```

Use mixed when iCMS2 and iCMS3 scores are close and malignant evidence is clear. Use indeterminate when UMI/stress/doublet/marker evidence is insufficient.

---

## 6. Functional overlays

| overlay | markers | final suffix | notes |
|---|---|---|---|
| Stem-like / crypt-base | `LGR5, OLFM4, ASCL2, SMOC2, AXIN2, RNF43, LRIG1, PROM1, CD44, SOX9` | `StemLike` | requires malignant evidence |
| TA / proliferative | `MKI67, TOP2A, PCNA, MCM2-7, TYMS, CENPF, UBE2C, CDK1, CCNB1` | `TAProliferative` | cycling overlay; can dominate clustering |
| Enterocyte / absorptive-like | `CA1, CA2, AQP8, SLC26A3, CEACAM7, MS4A12, GUCA2A/B, KRT20, FABP1, APOA1/4` | `AbsorptiveLike` | easily confused with normal colonocyte |
| Goblet / mucinous | `MUC2, TFF3, SPINK4, CLCA1, FCGBP, AGR2, REG4, MUC5AC, MUC5B` | `GobletMucinous` | require multiple mucin genes; check ambient mucin RNA |
| BEST4+ absorptive-like | `BEST4, OTOP2, SPIB, GUCA2B, CA7, NPY` | `BEST4Like` | differentiation identity only; adjudicate malignancy using CNV/SNV, normal atlas and patient/tissue context |
| Enteroendocrine-like | `CHGA, CHGB, NEUROD1, PAX6, PYY, GCG, TPH1, SCG5` | `EnteroendocrineLike` | if `SYP/INSM1/DLL3/ASCL1`, consider neuroendocrine component |
| Tuft-like | `DCLK1, TRPM5, POU2F3, AVIL, SH2D6, GFI1B` | `TuftLike` | rare; exclude normal tuft |
| Paneth-like / DCS | `DEFA5, DEFA6, LYZ, REG3A, PLA2G2A, MMP7, KIT, SOX9` | `PanethLike` | `LYZ` alone insufficient |
| Inflammatory / IFN | `IFIT1/2/3, ISG15, MX1, OAS1/2, STAT1, IRF1, CXCL9/10/11, HLA-A/B/C, B2M, TAP1` | `InflammatoryIFN` | exclude immune doublet |
| Metabolic | `SLC2A1, LDHA, HK2, PGK1, FABP1, APOA1/4, MTTP, ALDOB` | `Metabolic` | use cNMF/pathway context |
| EMT / invasive / tumor-budding | `LAMC2, LAMB3, ITGA5, ITGB1, MMP7, MMP14, PLAU, SERPINE1, TGFBI, EMP1, L1CAM, VIM, S100A4, ZEB1/2, SNAI2` | `EMTInvasive` | must retain epithelial/malignant evidence; exclude CAF |
| Fetal / regenerative / oncofetal | `CLU, ANXA1, L1CAM, TACSTD2, KRT17, KRT19, SOX9, EMP1, AREG, EREG, FOS, JUN` | `FetalRegenerative` | treatment/plasticity context helps |
| Hypoxia / stress / UPR | `CA9, VEGFA, SLC2A1, LDHA, ENO1, PGK1, NDRG1, BNIP3, ADM, HSPA1A/B, DNAJB1, ATF3, DDIT3, XBP1, HSPA5` | `HypoxiaStress` | if all lineages high, suspect dissociation artifact |
| Drug-tolerant / persister-like | `MEX3A, EMP1, L1CAM, CLU, ANXA1, KRT17, ALDH genes, ABCB1, ABCG2`; low `MKI67/LGR5-WNT` | `DTP-like` | requires treatment or strong context; otherwise provisional |

---

## 7. Orthogonal classifiers

```text
CMS1-4: sample / region pseudobulk only, not per-cell epithelial labels.
CRIS-A-E: clone / cluster / pseudobulk auxiliary scores, not sparse single-cell hard labels.
IMF: sample / region label combining iCMS + MSI + fibrosis.
PDS/SMI: useful for iCMS2/CMS2 internal structure; not a replacement for iCMS.
```

Do not write `CMS4 cell`. Correct phrasing:

```text
sample-level CMS4-like/fibrotic; malignant epithelial cells include EMT-invasive/pEMT state.
```

---

## 8. Final naming examples

```text
CRC_MalignantEpi_iCMS2_StemLike_CNVhigh
CRC_MalignantEpi_iCMS2_TAProliferative_CNVhigh
CRC_MalignantEpi_iCMS3_GobletMucinous_MSIlike
CRC_MalignantEpi_iCMS3_InflammatoryIFN_MSI
CRC_MalignantEpi_iCMS3_EMTInvasive_FibroticRegion
CRC_MalignantEpi_iCMSmixed_FetalRegenerative_DTP-like
CRC_NormalEpi_AbsorptiveColonocyte
CRC_NormalEpi_BEST4Absorptive
CRC_NormalEpi_Goblet
CRC_AmbiguousEpi_CNVnegative_TumorRegion
```

---

## 9. do_not_call

```text
Do not call CMS1/2/3/4 at single-cell epithelial level.
Do not call malignant stem/TA without CNV/SNV/tumor-cluster evidence.
Do not call iCMS3 CNV-neutral cells normal by default; MSI/iCMS3 can be CNV-low.
Do not call EMT tumor if epithelial and malignant evidence are absent; consider CAF.
Do not call goblet/mucinous from single MUC2 gene.
Do not call BEST4/absorptive/goblet malignant unless malignant evidence exists.
Do not call DTP-like without treatment or strong persister signature context.
```
