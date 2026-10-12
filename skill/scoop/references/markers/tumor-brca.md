# BRCA 肿瘤细胞注释 skill：乳腺癌恶性上皮与正常乳腺实质上皮

版本：v1.0-refined  
范围：human breast tumor scRNA-seq。  
定位：只在 `Epithelial/parenchymal` 候选细胞上使用；必须先经过 FastCNVpy / SNV / reference 判定 malignant_status。

---

## 1. 决策流程

```text
1. Epithelial gate: EPCAM/TACSTD2/KRT8/KRT18/KRT19/CDH1/MUC1
2. Exclude immune, endothelial, fibroblast, mural, adipocyte doublets
3. Assign malignant_status using FastCNVpy/SNV/reference/spatial/pathology
4. CNV-low/neutral cells -> joint malignancy review; normal breast phenotypes require normal-reference and tissue-context support
5. malignant cells -> BRCA malignant intrinsic subtype
6. cNMF programs -> functional overlays: HRhi, Prolifhi, HER2Signalhi, BasalStemhi, EMThi, Hypoxiahi, IFN_MHChi, etc.
```

---

## 2. Breast epithelial gate

| category | markers |
|---|---|
| epithelial positive | `EPCAM, TACSTD2, KRT8, KRT18, KRT19, CDH1, MUC1, CLDN3, CLDN4` |
| exclude immune | `PTPRC, CD3D, MS4A1, LST1, CD68` |
| exclude endothelial | `PECAM1, VWF, KDR, CLDN5` |
| exclude fibroblast/CAF | `COL1A1, DCN, LUM, PDGFRA, FAP, POSTN` |
| exclude mural | `RGS5, CSPG4, PDGFRB, MYH11, MCAM` |
| exclude adipocyte | `ADIPOQ, PLIN1, FABP4` |

---

## 3. 正常乳腺实质上皮细胞

Normal origin requires coherent breast lineage plus normal-reference and tissue-context support. CNV neutrality or absent/unassayed tumour SNVs alone does not establish normality.

| phenotype | markers | final label | notes |
|---|---|---|---|
| Luminal epithelial general | `EPCAM, KRT8, KRT18, KRT19, CDH1, MUC1` | `BRCA_NormalEpi_Luminal` | parent luminal label |
| Luminal hormone-sensing cell | `ESR1, PGR, AR, FOXA1, GATA3, TFF1, TFF3, BCL2` | `BRCA_NormalEpi_LHS` | low proliferation expected |
| Luminal adaptive secretory progenitor / luminal progenitor | `EHF, ELF5, KIT, ALDH1A3, LTF, SLPI, S100A14, KRT8, KRT18, KRT19` | `BRCA_NormalEpi_LASP` | `ESR1/PGR` low or negative |
| Mature secretory / lactocyte-like | `LALBA, CSN2, CSN3, LTF, WFDC2, ELF5, MUC1` | `BRCA_NormalEpi_Secretory` or `BRCA_NormalEpi_Lactocyte` | depends on pregnancy/lactation/background |
| Basal / myoepithelial cell | `KRT5, KRT14, KRT15, KRT17, TP63, ACTA2, MYL9, TAGLN, ITGA6` | `BRCA_NormalEpi_BasalMyo` | normal duct/lobule bilayer context supports |
| Cycling epithelial | parent marker + `MKI67, TOP2A, UBE2C, CDC20, CDK1, CCNB1, PCNA, MCM2, MCM5, MCM6` | `BRCA_NormalEpi_Cycling_<parent>` | normal-reference/context support required; cycling does not imply malignant |

---

## 4. BRCA malignant epithelial subtypes

Only use the following labels for `malignant-high` or `malignant-probable` epithelial cells.

| subtype | markers / signature | negative / caution | final label |
|---|---|---|---|
| Luminal A-like | `ESR1, PGR, FOXA1, GATA3, XBP1, TFF1, TFF3, BCL2, PIP, SCGB2A2, AGR2, AGR3`; low proliferation | `MKI67/TOP2A/UBE2C/CCNB1` low; basal low; HER2 low | `BRCA_MalignantEpi_LuminalA` |
| Luminal B-like | luminal markers plus high `MKI67, TOP2A, UBE2C, CDC20, CDK1, CCNB1, CCNB2, AURKB` | basal low; HER2 variable | `BRCA_MalignantEpi_LuminalB` |
| HER2-enriched | `ERBB2, GRB7, STARD3, MIEN1, PGAP3, TCAP, KRT7, KRT19`; 17q12 amplification if available | `ESR1/PGR` often low/intermediate but not always absent | `BRCA_MalignantEpi_HER2E` |
| Basal-like / TNBC-like | `KRT5, KRT14, KRT17, KRT6A, KRT6B, KRT6C, EGFR, TP63, ITGA6, LAMC2, S100A2, CD44` | `ESR1/PGR/ERBB2` low; clinical TNBC and basal-like are related but not identical | `BRCA_MalignantEpi_BasalLike` or `BRCA_MalignantEpi_TNBC_BasalLike` |
| Claudin-low / mesenchymal-like | high `VIM, FN1, ZEB1, ZEB2, SNAI2, TWIST1, AXL, CDH2, COL1A1, COL1A2, ACTA2, SPARC, ITGA5, CD44, ALDH1A1`; low `CLDN3, CLDN4, CLDN7, EPCAM, CDH1, MUC1, CD24` | best used as phenotype/tag, not mutually exclusive subtype; exclude CAF doublet | `BRCA_MalignantEpi_ClaudinLowLike` or overlay `EMThi` |
| TNBC-LAR / AR-luminal-like | `AR, FOXA1, GATA3, XBP1, KRT8, KRT18, MUC1`; usually low ER/PR/HER2 | use mainly in TNBC/HR-low malignant cells | `BRCA_MalignantEpi_TNBC_LAR` or `BRCA_MalignantEpi_BasalLike_ARhi` |
| Normal-like | not recommended as primary tumor label | often normal admixture, low purity, doublet, or ambient RNA | use `BRCA_AmbiguousEpi_NormalLike_CNVreview` pending joint genomic, reference and patient/tissue review |

---

## 5. cNMF / functional overlays

These overlays can be attached to any malignant subtype when supported by cNMF and marker evidence.

| overlay | markers | label suffix |
|---|---|---|
| Hormone response / ER signaling | `ESR1, PGR, FOXA1, GATA3, XBP1, TFF1, TFF3, BCL2, PIP, SCGB2A2` | `HRhi` |
| Proliferative / cycling | `MKI67, TOP2A, UBE2C, CDC20, CDK1, CCNB1, CCNB2, AURKB, MCM2, MCM5, MCM6, PCNA` | `Prolifhi` |
| HER2 / RTK activated | `ERBB2, GRB7, STARD3, PGAP3, MIEN1, EGFR, ERBB3, AREG` | `HER2Signalhi` |
| Basal-stem / squamous-like | `KRT5, KRT14, KRT17, KRT6A/B/C, TP63, EGFR, ITGA6, LAMC2, S100A2, CD44` | `BasalStemhi` |
| EMT / claudin-low | high `VIM, FN1, ZEB1/2, SNAI2, TWIST1, AXL, CDH2, COL1A1/2, ACTA2, SPARC`; low epithelial junction genes | `EMThi` or `ClaudinLowhi` |
| Hybrid epithelial-mesenchymal plasticity | simultaneous `EPCAM/CDH1/KRT8/18` and `VIM/FN1/S100A2/CRYAB/KRT14` | `HybridEMP` |
| Hypoxia / glycolysis | `CA9, VEGFA, HILPDA, NDRG1, BNIP3, LDHA, SLC2A1, ENO1, PGK1` | `Hypoxiahi` |
| Stress / UPR | `HSPA1A, HSPA1B, HSP90AA1, HSPA5, DDIT3, ATF4, XBP1, DNAJB1, HSPB1` | `StressUPRhi` |
| Antigen presentation / IFN | `HLA-A/B/C, B2M, TAP1/2, NLRC5, HLA-DRA/DRB1/DPA1/DPB1, CD74, STAT1, IRF1, ISG15, IFI6, IFI27, IFIT1/3, CXCL9/10/11` | `IFN_MHChi` |
| Secretory / lactation-like | `PIP, SCGB2A2, LTF, SLPI, WFDC2, ELF5, AGR2, AGR3, MUC1` | `Secretoryhi` |
| AR / LAR | `AR, FOXA1, GATA3, XBP1, KRT8, KRT18, MUC1` | `ARhi` or `LAR` |
| Inflammatory cytokine / chemokine | `CXCL1, CXCL2, CXCL3, CXCL8, IL6, CCL2, CSF1, CXCL16` | `InflamChemohi` |

---

## 6. Final naming examples

```text
BRCA_NormalEpi_LHS
BRCA_NormalEpi_LASP
BRCA_NormalEpi_BasalMyo
BRCA_NormalEpi_Cycling_LASP
BRCA_MalignantEpi_LuminalA_HRhi_CNVhigh
BRCA_MalignantEpi_LuminalB_Prolifhi_CNVhigh
BRCA_MalignantEpi_HER2E_HER2Signalhi_CNVhigh
BRCA_MalignantEpi_BasalLike_BasalStemhi_CNVhigh
BRCA_MalignantEpi_BasalLike_EMThi_CNVhigh
BRCA_MalignantEpi_TNBC_LAR_ARhi_CNVhigh
BRCA_AmbiguousEpi_CNVweak
```

---

## 7. do_not_call

```text
Do not call malignant breast epithelial from EPCAM/KRT alone.
Do not call Luminal A vs Luminal B without proliferation evidence.
Do not call HER2E from ERBB2 single gene without GRB7/STARD3/amp or module support.
Do not call claudin-low as a mutually exclusive subtype unless the whole signature supports it.
Do not assign malignancy from normal-like expression alone; review normal admixture, ambient RNA and doublets alongside genomic and patient/tissue evidence. CNV-low alone does not exclude malignancy.
Do not interpret cNMF state overlays as primary tumor subtype.
```
