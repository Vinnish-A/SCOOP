# T 细胞注释 skill：历史表型 + 2025 modular nomenclature

版本：v1.0-refined  
范围：human / mouse T cells；优先适配 human tumor scRNA-seq。  
信息来源：T 细胞命名遵循 Masopust 等 2025 年 T cell nomenclature consensus statement 的原则：明确实验依据、使用既有术语字典、并在需要时采用 modular nomenclature 表示多个独立属性。  
本文件只提供真实可追溯的 T 细胞表型与 marker；功能连续谱由 cNMF / FastCNMF 提供。

---

## 1. 决策流程

```text
1. T gate: CD3D/CD3E/CD3G/TRAC/TRBC/CD2/CD7
2. Exclude pure NK, B, myeloid, epithelial, stromal doublet
3. Determine lineage: CD4 / CD8 / gamma-delta / MAIT / iNKT / NKT-like / unresolved T
4. Determine historical phenotype: naive, activated, memory, effector, cytotoxic, helper, Treg, TFH, exhausted, TRM, etc.
5. Determine migration and differentiation attributes only when markers/context support them
6. Attach cNMF-derived functional overlays: cytotoxic, IFN, cycling, ribosomal, hypoxia, stress, etc.
7. For uncertain cells, use real phenotype + function rather than inventing a new subtype
```

---

## 2. T cell gate

| item | markers |
|---|---|
| transcript positive | `CD3D, CD3E, CD3G, TRAC, TRBC1, TRBC2, CD2, CD7` |
| protein positive | `CD3, TCRαβ, TCRγδ, CD2, CD7` |
| minimum call | `CD3D/CD3E/TRAC` positive or protein `CD3` positive |
| exclude NK without TCR | `KLRD1, KLRF1, NCAM1, FCGR3A, NKG7, GNLY` with `TRAC/CD3D` negative |
| exclude B | `MS4A1, CD79A, CD79B` |
| exclude myeloid | `LYZ, LST1, FCN1, C1QA` |

---

## 3. Modular T 输出语法

```text
<lineage> T<function><migration><migration_subscript><differentiation><differentiation_subscript><antigen>
```

| field | allowed values | evidence |
|---|---|---|
| lineage | `CD4+`, `CD8+`, `gamma_delta+`, `MAIT`, `iNKT`, `NKT` | co-receptor/TCR markers |
| function | `TH1`, `TH2`, `TH9`, `TH17`, `TFH`, `Treg`, `TCTL`, `TC17`, `unknown` | function markers |
| migration | `S`, `D`, `U` | `S`: SLO-homing, `CCR7/SELL` high; `D`: disseminated, `CCR7/SELL` low; `U`: unknown |
| migration subscript | `B`, `W`, `R` | `B`: blood only; `W`: widespread recirculation; `R`: resident, requires tissue context |
| differentiation | `N`, `A`, `M`, `X`, `G` | naive, activated, memory, exhausted, anergic |
| differentiation subscript | `p`, `t` | progenitor / terminal |
| antigen | `0`, `+`, omit | antigen cleared/irrelevant, persistent antigen, unknown |

Recommended canonical project labels can be simpler:

```text
T_CD4_Naive
T_CD8_TCM
T_CD8_TRM
T_CD8_TXp
T_Treg_eTreg
T_MAIT
```

The modular code should be stored as an additional field, not necessarily used as the plot label.

---

## 4. Differentiation state panel

| phenotype | human markers / assays | mouse markers / assays | final label |
|---|---|---|---|
| Naive T cell | `CCR7hi, CD45RA+, CD45RO-, CD95-`; other `CD27+, CD28+, IL7R+, CD11Alo, CD122-`; poor short-term `IFNG/IL4/IL5` | `CD62L+, CD44low, CD11alow, CD122low`; poor short-term `IFNG/IL4/IL5` | `T_<lineage>_Naive` or modular `TN` |
| Activated T cell | transient `CD69+, IL2RA/CD25hi, KLF2low`; often `MKI67+, HLA-DRA/DRB1+, CD38hi` | transient `CD69+, CD25+, KLF2low`; often `MKI67+` | `T_<lineage>_Activated` or `TA` |
| Memory T cell | difficult to separate from other antigen-experienced T; `CD11Ahi`, often `CD95+, CD58+, CD49D+` | `CD44hi, CD11Ahi`, typically `IL7R/CD127+` | `T_<lineage>_Memory` or `TM` |
| Anergic T cell | poorly characterized | `CD44+, FR4+, NT5E/CD73+, FOXP3-, NRP1/CD304+, CXCR5-` | `T_Anergic` or `TG`; human low confidence |
| Virtual memory T cell | difficult to define in humans | `CD44hi, CD49dlow, CD62Lhi, CD122hi` | `T_VirtualMemory`; mouse only unless strong evidence |

---

## 5. CD4 functional phenotypes

| phenotype | positive markers | required negative / caution | final label |
|---|---|---|---|
| TH1 | `TBX21, CXCR3, IFNG` | not `IL4/IL5/IL13` dominated | `T_CD4_TH1` |
| TH2 | `GATA3, PTGDR2/CRTH2, CCR4, IL4, IL5, IL13` | not `IFNG` dominated | `T_CD4_TH2` |
| TH9 | `SPI1/PU.1, IRF4, BATF, CCR8, IL9` | rare in scRNA; require multi-marker support | `T_CD4_TH9` |
| TH17 | `RORC, CCR6, IL17A, IL17F` | distinguish from MAIT/γδ IL17 programs | `T_CD4_TH17` |
| CD4 CTL | `CD4, GZMB, PRF1, NKG7, GNLY`, often `CCR7-, CD27-` | exclude NK/T doublet; need CD3/TRAC | `T_CD4_CTL` |
| TFH | `CXCR5, BCL6, ICOS, PDCD1` | GC/TFR localization requires tissue/imaging context | `T_CD4_TFH` |
| GC TFH | `BCL6hi, CXCR5hi, PDCD1hi, ICOS+`; mouse also `FR4hi, PSGL1low` | use only with lymphoid/follicle context | `T_CD4_GC_TFH` |
| TFH memory | human `CXCR4+, FOXP3-`; mouse `TCF7+, FR4+, P2RX7+, CXCR4+` | after antigen clearance only if supported | `T_CD4_TFHM` |
| Treg | `FOXP3, IL2RA/CD25, CTLA4, TNFRSF18/GITR`, `IL7Rlow` | `FOXP3` alone insufficient | `T_Treg` |
| eTreg | `FOXP3hi, CCR8hi, ICOS, IRF4, IL10, TGFB1`; human often `CD45RO+` | state overlay of Treg | `T_Treg_eTreg` |
| tTreg / pTreg | human hard to distinguish; mouse tTreg `NRP1hi, IKZF2/Helioshi`, pTreg often `NRP1low, Helioslow/int` | human should usually remain `T_Treg` | `T_Treg_tTreg` or `T_Treg_pTreg` only when evidence supports |
| TFR | `CXCR5, FOXP3, GITR`, mouse also `PDCD1, BCL6` | must have follicle/GC localization support | `T_Treg_TFR` |

---

## 6. CD8 effector / cytotoxic / exhausted phenotypes

| phenotype | markers | final label | notes |
|---|---|---|---|
| Cytotoxic CD8 T | `CD8A, CD8B, NKG7, GZMB, PRF1, GNLY, IFNG` | `T_CD8_Cytotoxic` | cytotoxic is function/state; can coexist with memory/exhaustion |
| SLEC / terminal effector | human `CCR7-, CD45RA+, KLRG1+, TCF7low, TBX21hi, CX3CR1+, IL7Rlow, CD27-, B3GAT1/CD57+`; mouse `KLRG1+, IL7Rlow` | `T_CD8_SLEC` or `T_CD8_TE` | no definitive marker; call cautiously |
| MPEC | mouse `IL7R/CD127+, CD27+, KLRG1low, TCF7+, TBX21int, CD62L+/-, CX3CR1-, CXCR3+`; human markers not universal | `T_CD8_MPEC` | use mainly in mouse or well-controlled datasets |
| EEC | `KLRG1low, IL7Rlow, CXCR3, TBX21, TCF7low, IFNG, TNF` | `T_CD8_EEC` | activated effector before SLEC/MPEC split |
| LLEC | `KLRG1hi, IL7Rint, CX3CR1hi, GZMBhi` | `T_CD8_LLEC` | long-lived cytotoxic memory-like effector |
| TPEX / progenitor exhausted | `PDCD1, TOX, TCF7, BCL6, SLAMF6, CXCR3, LEF1, CD28, NT5E, XCL1, CXCR5`; low/negative `HAVCR2, ENTPD1, CX3CR1hi, GZMBhi` | `T_CD8_TXp` | PDCD1 alone is not enough |
| Intermediate / effector exhausted | `PDCD1, TOX, HAVCR2/TIM3, TBX21, GZMB, PRF1, IFNG, CX3CR1`; low `TCF7, SLAMF6, CD101` | `T_CD8_TX` or `T_CD8_TXeff` | often cytotoxic + exhausted program |
| Terminal exhausted | `PDCD1, TOX, HAVCR2/TIM3, GZMB, ENTPD1/CD39, CD244/2B4, CD101`; low `TCF7, SLAMF6, CX3CR1, CXCR3` | `T_CD8_TXt` | persistent antigen context increases confidence |

---

## 7. Memory / migration phenotypes

| phenotype | human markers | mouse markers | final label |
|---|---|---|---|
| TCM | `SELL/CD62L+, CCR7+, CD27+, CX3CR1-, CD45RO+, CD45RA-` | `CD44hi, CD62L+`, often `CD69-, CD11ahi, IL7R+` | `T_<lineage>_TCM` |
| TEM | `CCR7-, CD45RA-`, often `SELL/CD62L-` | `CD44hi, CD62L-` | `T_<lineage>_TEM` |
| TEMRA | human only: `CCR7-, CD45RA+, CD27-, CX3CR1hi, B3GAT1/CD57+`, cytolytic | not applicable | `T_CD8_TEMRA` |
| TRM | `KLF2low, CD69+, SELL/CD62L-, CCR7-`; often `ITGAE/CD103+, ITGA1/CD49a+, CXCR6+, CD101+` | similar; often `P2RX7+` | `T_<lineage>_TRM`; requires tissue context |
| Ex-TRM | circulating cell with prior TRM evidence | mouse may retain `CD43/1B11` | `T_<lineage>_exTRM`; rare, needs direct evidence |
| TPM | `CD45RA-, CCR7int, CD27+` | `CX3CR1int, CCR7int, CD27+` | `T_<lineage>_TPM`; no definitive marker |
| TSCM | `CD95+, CCR7+, CD27+, CD28+, CD45RA+, EOMES-` | `CD44low, CD62L+, CD122+, CD27+, Sca1+, IL7R+` | `T_<lineage>_TSCM` |

---

## 8. Innate-like T phenotypes

| phenotype | markers / gold standard | final label | notes |
|---|---|---|---|
| MAIT | `CD3+, TRAC+, KLRB1/CD161+, SLC4A10, ZBTB16, TRAV1-2`; gold standard `MR1-5OPRU tetramer+` | `T_MAIT` | blood flow proxy: `CD3+ CD161+ Vα7.2+` |
| iNKT | `CD3+, TRAC+, ZBTB16`, `TRAV10-TRAJ18` proxy; gold standard `CD1d-αGalCer tetramer+` | `T_iNKT` | blood proxy: `Vα24Jα18` or `Vα24/Vβ11` |
| γδ T | `TRDC, TRGC1, TRGC2`; protein `γδTCR+` | `T_gamma_delta` | human blood Vγ9Vδ2 proxy: `TRGV9, TRDV2, KLRB1, KLRC1` |

---

## 9. cNMF overlays for T cells

The following labels are overlays, not standalone T cell types:

```text
Cytotoxic
IFNhi
Cycling
RibosomalHigh
Stress
TRM_program_high
Exhaustion_program_high
APC_like_MHCIIhi
```

Example final records:

```yaml
- final_label: T_CD8_TXp
  modular_code: CD8+ TUXp+
  functional_overlay: none
- final_label: T_CD8_Cytotoxic
  modular_code: CD8+ TCTLDM
  functional_overlay: Cycling
- final_label: T_CD4_Treg_eTreg
  modular_code: CD4+ TregUDA
  functional_overlay: IFNhi
```

---

## 10. do_not_call

```text
Do not call exhausted T from PDCD1 alone.
Do not call Treg from FOXP3 alone.
Do not call TRM without tissue context or residence-associated evidence.
Do not call TFR without follicle/GC localization support.
Do not call TEMRA in mouse.
Do not force TH1/TH2/TH17 labels when cytokine genes are absent and only chemokines are weakly present.
Do not create new T cell types for cycling, ribosomal, stress, or hypoxia programs.
```
