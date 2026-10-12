# 通用细胞命名策略：属 + 种差 + 状态后缀

版本：v1.0-refined

---

## 1. 核心思想

推荐采用“**属 + 种差**”结构：

```text
<Context>_<Genus>_<HistoricalPhenotype>_<SpecificDifference>_<OptionalState/Evidence>
```

- `Genus`：最稳定的细胞属类，例如 `T`、`B`、`Mac`、`CAF`、`Endo`、`MalignantEpi`。
- `HistoricalPhenotype`：文献中已存在、可由 marker 支持的真实表型，例如 `CD8_TEMRA`、`Treg`、`SPP1_TREM2_TAM`、`iCMS2`、`AC-like`。
- `SpecificDifference`：同属内最关键的种差，例如 `LuminalA` vs `BasalLike`、`iCMS2` vs `iCMS3`、`TRM-like` vs `Mo-like`。
- `State/Evidence`：由 cNMF、CNV、空间或病理支持的附加说明，例如 `IFNhi`、`Hypoxiahi`、`CNVhigh`、`Cycling`。

---

## 2. canonical label 与 display label

建议同时保留长标签和短标签：

```yaml
canonical_label: CRC_MalignantEpi_iCMS3_GobletMucinous_CNVevident
short_label: iCMS3 Goblet
plot_label: CRC iCMS3 goblet
```

`canonical_label` 用于表格、复核和跨项目合并；`short_label` 用于 UMAP、dotplot、heatmap。

---

## 3. 推荐命名语法

### 3.1 免疫细胞

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

### 3.2 基质细胞

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

### 3.3 肿瘤细胞

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

### 3.4 正常实质细胞

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

## 4. 状态后缀白名单

以下后缀只能作为状态或 overlay，不能单独作为 cell type：

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

## 5. 合并规则

两个小群应合并同名的情况：

```text
same broad lineage
same historical phenotype
same key exclusion markers
same malignant/normal status if tumor
only different in cycling/ribosomal/stress/mitochondrial/cNMF generic program
```

示例：

```text
T_CD8_Cytotoxic + T_CD8_Cytotoxic_Cycling -> figures 可分开显示，但 final parent label 同为 T_CD8_Cytotoxic
Mac_TAM_SPP1_TREM2 + Mac_TAM_SPP1_TREM2_Cycling -> parent label 同为 Mac_TAM_SPP1_TREM2
CRC_MalignantEpi_iCMS2_TAProliferative + CRC_MalignantEpi_iCMS2_TAProliferative_RibosomalHigh -> parent label 同为 CRC_MalignantEpi_iCMS2_TAProliferative
```

---

## 6. 不确定命名

### 6.1 真实身份明确，精细表型不明确

```text
T_CD4_Unresolved
Mac_Unresolved
CAF_Unresolved
MalignantEpi_Unresolved
```

### 6.2 两个历史表型都可能成立

```text
T_CD8_TEXp_or_TCM_lowConfidence
Mac_C1QC_or_SPP1_intermediate
CRC_MalignantEpi_iCMS2_iCMS3_mixed
GBM_Malignant_AC_MES_hybrid
```

### 6.3 肿瘤/正常状态不确定

```text
BRCA_AmbiguousEpi_CNVweak
CRC_AmbiguousEpi_CNVnegative_TumorRegion
GBM_PutativeMalignant_OPC-like_lowConfidence
```

---

## 7. 禁止命名

```text
禁止：T_Pdcd1_exhausted                 # PDCD1 单基因不足
禁止：Mac_M1 / Mac_M2                    # M1/M2 不作为最终标签
禁止：CAF_FAP                            # FAP 单基因不足
禁止：CMS4_cell                          # CMS 是 bulk/sample/TME 层，不是单细胞上皮标签
禁止：GBM_MES_only_by_CD44               # CD44 单基因不足
禁止：BRCA_NormalLike_malignant          # normal-like 本身不证明恶性；联合基因组、正常图谱、患者/组织证据复核，不要求CNV必须升高
禁止：Cycling_cell                        # cycling 是状态，不是独立身份
```
