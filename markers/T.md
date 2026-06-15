# 快速、稳健、简约的 T 细胞注释

版本：v0.1  
原则：只声明已被数据支持的属性；用正标记 + 排除标记 + 组织背景 + CNV/参考图谱交叉验证；证据不足时输出 `others*`，不要强行细分。

---

## 0. Agent 通用执行规则

### 0.1 输入
- 表达矩阵：gene × cell / protein × cell。
- 聚类结果：cluster id、DE genes、marker scores。
- 可选：CNV 推断、TCR/BCR、CITE-seq protein、空间坐标、样本临床信息。

### 0.2 输出 schema
```yaml
cell_annotation:
  broad_type: <major compartment>
  subtype: <cell subtype, optional>
  state: <activation/differentiation/stress/cycling/etc, optional>
  modular_code: <for T cells or other modular labels, optional>
  positive_evidence: [genes/proteins/assays]
  exclusion_evidence: [negative markers or absent compartments]
  context_evidence: [CNV, tissue, spatial, reference atlas]
  confidence: high | medium | low
  uncertainty: <why ambiguous or which alternatives remain>
  do_not_call: [labels that should not be used with current evidence]
```

### 0.3 置信度
- `high`：≥2 个核心正标记 + 关键排除标记成立 + 与组织背景/参考图谱一致。
- `medium`：1 个核心标记强阳性，或多标记弱阳性，但排除项不足。
- `low`：marker 稀疏、污染/ambient RNA 可能性高、或相邻谱系难以区分。

### 0.4 禁止规则
- 不凭单个 marker 注释精细亚群。例如 `FOXP3` 单独不足以稳健判定 Treg；`PDCD1` 单独不足以判定 exhausted T。
- 不把所有 `CD44+` 肿瘤细胞都叫 MES-like；必须结合 MES 程序或多基因 score。
- 不把所有 `P2RY12/TMEM119` 下降的 myeloid 都叫 macrophage；GBM 微胶质细胞会失去部分稳态 marker。
- 不用 M1/M2 作为最终 macrophage 标签；最多作为外部比较分数。
- 不把 `HBB/HBA/PPBP/PF4` 高表达细胞作为真正组织细胞，除非有空间/蛋白证据。

---

# 1. T 细胞注释 Skill

## 1.1 Skill ID
```yaml
skill_id: T_CELL_MODULAR_ANNOTATION_V1
scope: human and mouse T cells; optimized for human tumor tissue single-cell data
principle: modular; measured properties only; no overclaiming
```

## 1.2 T cell gate
```yaml
positive:
  transcript: [CD3D, CD3E, CD3G, TRAC, TRBC1, TRBC2, CD2, CD7]
  protein: [CD3, TCRab, TCRgd, CD2, CD7]
exclude:
  NK_without_TCR: [KLRD1, KLRF1, NCAM1, FCGR3A, NKG7, GNLY] with TRAC/CD3D negative
  B_cell: [MS4A1, CD79A, CD79B]
  myeloid: [LYZ, LST1, FCN1, C1QA]
minimum_call: CD3D/CD3E/TRAC positive or protein CD3 positive
```

## 1.3 T cell modular output grammar
```yaml
syntax: <lineage> T<function><migration><migration_subscript><differentiation><differentiation_subscript><antigen>
examples:
  - CD4+ TN
  - CD8+ TDRXp+
  - CD4+ TH17DM
  - CD8+ TCTLDA
fields:
  lineage: [CD4+, CD8+, gamma_delta+, MAIT, iNKT, NKT]
  function: [TH1, TH2, TH9, TH17, TFH, Treg, TCTL, TC17, TC22, unknown]
  migration:
    S: SLO-homing; CCR7/SELL high
    D: disseminated; CCR7/SELL low or absent
    U: unknown migration
  migration_subscript:
    B: blood-derived only; no tissue residence claim
    W: widespread non-lymphoid recirculation; require direct evidence
    R: resident; require tissue context plus CD69/ITGAE/CXCR6 or stronger evidence
  differentiation:
    N: naive
    A: activated
    M: memory
    X: exhausted
    G: anergic
  differentiation_subscript:
    p: progenitor/stem-like
    t: terminal
  antigen:
    '+': persistent antigen reasonably supported
    '0': antigen cleared or irrelevant
    '?': unknown; omit in final label if not assessed
```

## 1.4 Differentiation state panel
```yaml
naive_T:
  human_core: [CCR7, SELL, IL7R, TCF7, LEF1]
  protein_human: [CCR7hi, CD45RA+, CD45RO-, CD95-]
  mouse_protein: [CD62L+, CD44low, CD11alow, CD122low]
  output: TN
activated_T:
  core: [CD69, IL2RA, MKI67, CD38, HLA-DRA, HLA-DRB1]
  protein_human: [CD69+, CD25hi, KLF2low, Ki67+, HLA-DR+, CD38hi]
  output: TA
memory_T:
  core: [IL7R, CD44, CD69_or_CCR7_context_dependent]
  protein_human: [CD11ahi, CD95+, CD58+, CD49d+]
  output: TM
anergic_T:
  core_mouse: [FR4, NT5E, NRP1, CXCR5_negative, FOXP3_negative]
  caution: human markers poorly characterized
  output: TG
```

## 1.5 CD4 functional panel
```yaml
TH1:
  positive: [TBX21, CXCR3, IFNG]
  negative: [IL4, IL5, IL13]
TH2:
  positive: [GATA3, CCR4, PTGDR2, IL4, IL5, IL13]
  negative: [IFNG]
TH9:
  positive: [SPI1, IRF4, BATF, CCR8, IL9]
TH17:
  positive: [RORC, CCR6, IL17A, IL17F]
TFH:
  positive: [CXCR5, BCL6, ICOS, PDCD1]
  note: GC/TFR requires tissue localization or imaging support
Treg:
  positive: [FOXP3, IL2RA, CTLA4, TNFRSF18, IKZF2]
  negative_or_low: [IL7R]
  eTreg_positive: [FOXP3, CCR8, ICOS, IRF4, IL10, TGFB1]
CD4_CTL:
  positive: [CD4, GZMB, PRF1, NKG7, GNLY]
```

## 1.6 CD8 effector and exhausted panel
```yaml
cytotoxic_CD8:
  positive: [CD8A, CD8B, NKG7, GZMB, PRF1, GNLY, IFNG]
SLEC_or_terminal_effector:
  protein_human: [CCR7-, CD45RA+, KLRG1+, TCF1low, TBX21hi, CX3CR1+, IL7Rlow, CD27-, CD57+]
  transcript_proxy: [CX3CR1, KLRG1, GZMB, PRF1, TBX21, TCF7_low, IL7R_low]
MPEC_or_activated_progenitor:
  mouse_core: [IL7R, CD27, TCF7, CXCR3, KLRG1_low, CX3CR1_low]
EEC:
  positive: [KLRG1_low, IL7R_low, CXCR3, TBX21, TCF7_low, IFNG, TNF]
LLEC:
  positive: [KLRG1_hi, IL7R_int, CX3CR1_hi, GZMB_hi]
TPEX:
  positive: [PDCD1, TOX, TCF7, BCL6, SLAMF6, CXCR3, LEF1, CD28, NT5E, XCL1, CXCR5]
  negative: [HAVCR2, ENTPD1, CX3CR1_high, GZMB_high]
  output: TXp
TEX_int_or_TEX_eff:
  positive: [PDCD1, TOX, HAVCR2, TBX21, GZMB, PRF1, IFNG, CX3CR1]
  negative: [TCF7, SLAMF6, CD101]
  output: TX
TEX_term:
  positive: [PDCD1, TOX, HAVCR2, GZMB, ENTPD1, CD244, CD101]
  negative: [TCF7, SLAMF6, CX3CR1, CXCR3]
  output: TXt
```

## 1.7 Memory / migration panel
```yaml
TCM:
  human: [SELL, CCR7, CD27, CX3CR1_low, CD45RO, CD45RA_negative]
  modular: TSM
TEM:
  human: [CCR7_negative, CD45RA_negative, SELL_low]
  modular: TDM or TDBM if blood only
TEMRA:
  human_only: [CCR7_negative, CD45RA_positive, CD27_negative, CX3CR1_hi, B3GAT1/CD57]
TRM:
  positive: [CD69, ITGAE, ITGA1, CXCR6, KLF2_low, SELL_low, CCR7_low]
  require: tissue context; preferably spatial or prior residence evidence
  modular: TDRM
TSCM:
  human: [CD95, CCR7, CD27, CD28, CD45RA, EOMES_low]
  modular: TMp
TPM:
  human: [CD45RA_negative, CCR7_intermediate, CD27]
  caution: no definitive marker for non-lymphoid recirculation
```

## 1.8 Innate-like T panel
```yaml
MAIT:
  positive: [CD3D, TRAC, KLRB1, SLC4A10, ZBTB16, TRAV1-2]
  gold_standard: MR1-5OPRU tetramer+
iNKT:
  positive: [CD3D, TRAC, ZBTB16, TRAV10_TRAJ18_or_Va24-Ja18_proxy]
  gold_standard: CD1d-alphaGalCer tetramer+
gamma_delta_T:
  positive: [TRDC, TRGC1, TRGC2]
  human_blood_Vg9Vd2_proxy: [TRGV9, TRDV2, KLRB1, KLRC1]
```
