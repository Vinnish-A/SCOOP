# General Cell Naming Strategy: Genus + Differentia + State Suffix

Version: v1.0-refined

---

## 1. Core idea

Use the “**genus + differentia**” structure:

```text
<Context>_<Genus>_<HistoricalPhenotype>_<SpecificDifference>_<OptionalState/Evidence>
```

- `Genus`: the most stable cell genus, such as `T`, `B`, `Mac`, `CAF`, `Endo`, or `MalignantEpi`.
- `HistoricalPhenotype`: a real phenotype established in the literature and supportable by markers, such as `CD8_TEMRA`, `Treg`, `SPP1_TREM2_TAM`, `iCMS2`, or `AC-like`.
- `SpecificDifference`: the key differentia within the same genus, such as `LuminalA` vs `BasalLike`, `iCMS2` vs `iCMS3`, or `TRM-like` vs `Mo-like`.
- `State/Evidence`: supplementary information supported by cNMF, CNV, spatial evidence, or pathology, such as `IFNhi`, `Hypoxiahi`, `CNVhigh`, or `Cycling`.

---

## 2. Canonical label and display label

Retain both long and short labels:

```yaml
canonical_label: CRC_MalignantEpi_iCMS3_GobletMucinous_CNVevident
short_label: iCMS3 Goblet
plot_label: CRC iCMS3 goblet
```

Use `canonical_label` for tables, review, and cross-project merging; use `short_label` for UMAPs, dotplots, and heatmaps.

---

## 3. Recommended naming syntax

### 3.1 Immune cells

```text
T_CD4_Naive
T_CD8_TEMRA_Cytotoxic
T_CD8_TXp
T_CD8_TXt
T_Treg_eTreg
NK_Cytotoxic_FCGR3A
B_Naive
B_Memory
Plasma_IgG
Mac_Mo_FCN1_IL1B
Mac_TRM_FOLR2_LYVE1
Mac_TAM_SPP1_TREM2
DC_cDC1
DC_cDC2
DC_LAMP3_Migratory
```

### 3.2 Stromal cells

```text
FB_Resident_PI16
FB_Resident_DPT
CAF_Matrix_POSTN_COL11A1
CAF_myCAF_Contractile
CAF_iCAF_CXCL12_IL6
CAF_apCAF_MHCII
CAF_TLS_Reticular
Endo_Capillary
Endo_Venous
Endo_Lymphatic
Pericyte_RGS5
SMC_Contractile
```

### 3.3 Tumor cells

```text
BRCA_MalignantEpi_LuminalA_HRhi_CNVhigh
BRCA_MalignantEpi_LuminalB_Prolifhi_CNVhigh
BRCA_MalignantEpi_HER2E_HER2Signalhi_CNVhigh
BRCA_MalignantEpi_BasalLike_EMThi_CNVhigh
CRC_MalignantEpi_iCMS2_TAProliferative_CNVhigh
CRC_MalignantEpi_iCMS3_GobletMucinous_MSIlike
GBM_Malignant_AC-like_Gliosis
GBM_Malignant_OPC-like_Cycling
GBM_Malignant_MES-like_Hypoxia
```

### 3.4 Normal parenchymal cells

```text
BRCA_NormalEpi_LHS
BRCA_NormalEpi_LASP
BRCA_NormalEpi_BasalMyo
CRC_NormalEpi_AbsorptiveColonocyte
CRC_NormalEpi_Goblet
GBM_Normal_Astrocyte
GBM_Normal_OPC
GBM_Normal_Oligodendrocyte
```

---

## 4. State-suffix allowlist

The following suffixes may be used only as states or overlays, not as cell types on their own:

```text
Cycling
RibosomalHigh
Stress
HeatShock
IFNhi
APChi
Hypoxiahi
Glycolysishi
Prolifhi
Cytotoxic
EMThi
HybridEMP
DTP-like
Secretoryhi
HRhi
ARhi
HER2Signalhi
SPP1_TREM2
MHCIIhi
```

---

## 5. Merging rules

Two small groups should be merged under one name when:

```text
same broad lineage
same historical phenotype
same key exclusion markers
same malignant/normal status if tumor
only different in cycling/ribosomal/stress/mitochondrial/cNMF generic program
```

Examples:

```text
T_CD8_Cytotoxic + T_CD8_Cytotoxic_Cycling -> figures may show them separately, but the final parent label is the same: T_CD8_Cytotoxic
Mac_TAM_SPP1_TREM2 + Mac_TAM_SPP1_TREM2_Cycling -> the parent label is the same: Mac_TAM_SPP1_TREM2
CRC_MalignantEpi_iCMS2_TAProliferative + CRC_MalignantEpi_iCMS2_TAProliferative_RibosomalHigh -> the parent label is the same: CRC_MalignantEpi_iCMS2_TAProliferative
```

---

## 6. Uncertain naming

### 6.1 Real identity is clear, but the fine phenotype remains mixed after sufficient evidence is reviewed

If evidence has not yet been fully collected, record `evidence_pending` and actively obtain it; do not use the following unresolved names as final identities for missing evidence. This term is a work status, not a new runtime schema enum.

```text
T_CD4_Unresolved
Mac_Unresolved
CAF_Unresolved
MalignantEpi_Unresolved
```

### 6.2 Two historical phenotypes are both plausible

```text
T_CD8_TEXp_or_TCM_lowConfidence
Mac_C1QC_or_SPP1_intermediate
CRC_MalignantEpi_iCMS2_iCMS3_mixed
GBM_Malignant_AC_MES_hybrid
```

### 6.3 Tumor/normal status is uncertain

```text
BRCA_AmbiguousEpi_CNVweak
CRC_AmbiguousEpi_CNVnegative_TumorRegion
GBM_PutativeMalignant_OPC-like_lowConfidence
```

---

## 7. Prohibited names

```text
Do not use: T_Pdcd1_exhausted                 # PDCD1 alone is insufficient
Do not use: Mac_M1 / Mac_M2                    # M1/M2 are not final labels
Do not use: CAF_FAP                            # FAP alone is insufficient
Do not use: CMS4_cell                          # CMS is a bulk/sample/TME-level concept, not a single-cell epithelial label
Do not use: GBM_MES_only_by_CD44               # CD44 alone is insufficient
Do not use: BRCA_NormalLike_malignant          # normal-like alone does not establish malignancy; review jointly with genomic, normal-atlas, and patient/tissue evidence; elevated CNV is not required
Do not use: Cycling_cell                        # cycling is a state, not an independent identity
```
