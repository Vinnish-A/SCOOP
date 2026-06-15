下面是一套**可直接用于 CRC/COAD/READ 恶性上皮/腺癌细胞 scRNA-seq 注释**的分层策略。核心原则是：**不要把 bulk CMS 直接当作单细胞标签**；单细胞肿瘤上皮应以 **iCMS2/iCMS3 + 功能状态**为主轴，再叠加 CRIS、IMF、PDS/SMI、空间/多组学证据。经典 CMS 描述的是 bulk 肿瘤整体，CMS1–4 分别强调 MSI/immune、WNT-MYC/canonical、metabolic、mesenchymal/fibrotic/stromal；而 Nature Genetics 2022 的 CRC 单细胞/bulk 整合研究提出 iCMS2/iCMS3，并把 epithelial intrinsic state、MSI、fibrosis 解耦为 IMF 框架。([Nature][1])

---

## 1. 推荐输出标签体系

每个细胞最终输出 6 层标签：

```text
cell_id
├─ major_lineage: epithelial / immune / stromal / endothelial / unknown
├─ malignancy_status: malignant-high / malignant-probable / normal-epithelial / ambiguous
├─ tumor_intrinsic_subtype: iCMS2 / iCMS3 / iCMS-mixed / iCMS-indeterminate
├─ functional_state_primary: stem-like / TA-proliferative / absorptive / goblet-mucinous / inflammatory-IFN / EMT-invasive / fetal-regenerative / hypoxia-stress / DTP-persister / ...
├─ orthogonal_classifiers: CRIS-A/B/C/D/E score, PDS/SMI, CMS-pseudobulk, IMF sample/region class
└─ evidence/confidence: CNV, mutation, normal-reference mapping, module scores, spatial/pathology support
```

这样可以避免“CMS4 细胞”“CMS1 细胞”这类不严谨说法。更准确的写法应是：

```text
malignant epithelial cell; iCMS3; inflammatory-IFN functional state; sample-level CMS1-like; IMF=iCMS3_MSI
```

或：

```text
malignant epithelial cell; iCMS2; stem/TA-proliferative state; CRIS-D-like; sample-level CMS2-like
```

---

## 2. 命名来源与实际落地

| 命名系统                                                         |                                           来源层级 | 实际用途                                                                               | 不建议做法                                          |
| ------------------------------------------------------------ | ---------------------------------------------: | ---------------------------------------------------------------------------------- | ---------------------------------------------- |
| **CMS1–4**                                                   |                      bulk tumor / whole tissue | 对 patient/sample/region 做 pseudobulk 分类；解释免疫、WNT-MYC、metabolic、fibrotic/stromal 背景 | 不建议直接给单个细胞标 “CMS1 cell / CMS4 cell”            |
| **iCMS2 / iCMS3**                                            |   CRC 单细胞 + bulk 整合，tumor epithelial intrinsic | 恶性上皮细胞的主分类轴                                                                        | 不应忽略功能状态；iCMS2 内也有 differentiated/PDS3-like 细胞 |
| **IMF**                                                      |       epithelial intrinsic + MSI + fibrosis 解耦 | 样本/区域层面标签：iCMS2_MSS_NF、iCMS2_MSS_F、iCMS3_MSS_NF、iCMS3_MSS_F、iCMS3_MSI              | 不把 F 当作单个 EMT 细胞标签；F 需要 CAF/ECM/空间证据           |
| **CRIS-A–E**                                                 | CRC intrinsic subtype，尽量去除 stromal confounding | clone/pseudobulk/organoid/PDX 层面辅助解释                                               | 不建议对每个稀疏单细胞硬分 CRIS                             |
| **PDS / SMI**                                                |           2024 pathway-level subtype，尤其细分 CMS2 | 捕捉 iCMS2 内 canonical stem/proliferative 与 differentiated/slow-cycling 差异           | 不替代 iCMS，而是正交层                                 |
| **Epi-CRC / comparative atlas / perturbation-guided states** |                             2026 bioRxiv 预印本方向 | 作为扩展参考：hybrid/endoderm-like/oncofetal、扰动到机制映射                                      | 预印本未同行评议，不作为唯一金标准                              |

CRIS 五类可作为上皮内在辅助轴：CRIS-A 偏 mucinous/glycolytic/MSI/KRAS，CRIS-B 偏 TGFβ/EMT/poor prognosis，CRIS-C 偏 EGFR signaling，CRIS-D 偏 WNT/IGF2，CRIS-E 偏 Paneth-like/TP53。([Nature][2]) PDS/SMI 则把 PDS1 定义为 canonical/LGR5+ stem-rich/proliferative，PDS2 偏 regenerative/ANXA1+ stem-rich 并伴 stromal/immune TME，PDS3 是 CMS2 内 slow-cycling、differentiated lineage 增多且预后较差的一类；其分类器默认概率阈值为 0.6，低于阈值可标 mixed。([Nature][3])

---

## 3. 第一步：先识别上皮细胞，再判定恶性

### 3.1 上皮细胞筛选

**保留：**

```text
EPCAM, TACSTD2, KRT8, KRT18, KRT19, KRT20, CDH1, CLDN3, CLDN4, MSLN
```

**排除非上皮污染：**

```text
immune: PTPRC, CD3D, CD3E, TRAC, MS4A1, CD79A, LST1, FCGR3A, CD68
fibroblast/CAF: COL1A1, COL1A2, DCN, LUM, COL3A1, PDGFRA, ACTA2, TAGLN
endothelial: PECAM1, VWF, KDR, EMCN, PLVAP
pericyte/smooth muscle: RGS5, MCAM, CSPG4, MYH11
```

注意：EMT / tumor budding 细胞可能 **EPCAM 下降但 KRT8/KRT18/KRT19 仍阳性**。这类细胞不能因为 EPCAM 低就直接扔掉，但必须证明其有 epithelial + malignant 证据，否则很容易混入 CAF。

---

## 4. 正常实质上皮 vs 恶性上皮细胞

CRC 中正常实质细胞主要是正常肠上皮谱系。建议用 **CNV/突变/患者特异聚类/正常图谱映射/肿瘤程序**联合判定，而不是只靠 EPCAM 或 KRT。

### 4.1 正常肠上皮谱系 marker

| 正常谱系                               | marker                                                              |
| ---------------------------------- | ------------------------------------------------------------------- |
| Stem-like / crypt-base             | **LGR5, OLFM4, ASCL2, SMOC2, AXIN2, RNF43, LRIG1**                  |
| Transit-amplifying / proliferating | **MKI67, TOP2A, PCNA, MCM2, MCM5, TYMS, CENPF, UBE2C**              |
| Absorptive colonocyte              | **CA1, CA2, AQP8, SLC26A3, CEACAM7, MS4A12, GUCA2A, GUCA2B, KRT20** |
| Goblet                             | **MUC2, TFF3, SPINK4, CLCA1, FCGBP, AGR2**                          |
| BEST4+ absorptive                  | **BEST4, OTOP2, SPIB, GUCA2B, CA7, NPY**                            |
| Enteroendocrine                    | **CHGA, CHGB, NEUROD1, PAX6, PYY, GCG, TPH1, SCG5**                 |
| Tuft                               | **DCLK1, TRPM5, POU2F3, AVIL, SH2D6, GFI1B**                        |
| Paneth-like / deep crypt secretory | **DEFA5, DEFA6, LYZ, REG3A, PLA2G2A, MMP7, KIT**                    |

正常肠道 atlas 和健康结肠上皮单细胞图谱应作为 reference mapping 的基线；BEST4+、Paneth-like、stem/TA 等谱系在人体肠道中存在区域性和解剖部位差异，因此必须结合部位、邻近正常组织和病理背景解释。([HCA Data Explorer][4])

### 4.2 恶性判定证据

**A. CNV / aneuploidy 推断**

推荐工具：

```text
inferCNV / copyKAT / Numbat / HoneyBADGER
```

参考细胞：

```text
同一患者免疫细胞、内皮细胞、成纤维细胞、邻近正常上皮
```

CRC CIN 常见模式可作为支持证据：

```text
gains: 7p/q, 8q, 13q, 20p/q
losses: 1p, 4p/q, 8p, 14q, 15q, 17p, 18p/q
```

Nature Genetics 2022 的 iCMS/IMF 研究中，iCMS2 更常见这些染色体臂 CNV；iCMS3 可为 MSI-H 近二倍体，也可为 MSS 但 CNV 较少，因此 **CNV 阴性不能排除 iCMS3/MSI 型恶性细胞**。([Nature][1])

**B. 患者特异性突变**

优先使用 matched WES / panel / bulk RNA / targeted genotyping：

```text
APC, KRAS, NRAS, BRAF, TP53, PIK3CA, SMAD4, FBXW7, TCF7L2, ACVR2A, TGFBR2, MLH1/MSH2/MSH6/PMS2 status
```

scRNA 直接检突变灵敏度低，建议 cluster-level 汇总。只要某上皮 cluster 富集患者特异 driver mutation，即可强支持 malignant。

**C. 患者特异聚类**

恶性上皮细胞常呈患者特异 cluster；正常上皮细胞更容易跨患者按谱系混合。iCMS 研究中，正常上皮与肿瘤上皮的聚类行为不同，且该研究在 63 名患者中分析了大量上皮细胞并用 inferCNV 支持恶性识别。([Nature][5])

**D. 正常参考映射**

把上皮细胞映射到：

```text
matched adjacent normal
Human Cell Atlas gut / colon reference
Gut Cell Atlas
CRC atlas normal epithelial compartment
```

判定逻辑：

```text
normal epithelial:
  epithelial marker high
  CNV negative
  tumor SNV negative
  maps to normal lineage with high confidence
  shared across patients / adjacent normal

malignant-high:
  epithelial marker high
  CNV positive or tumor SNV positive
  patient-specific tumor cluster
  tumor epithelial program high

malignant-probable:
  epithelial marker high
  CNV negative or weak
  but tumor sample-specific cluster
  iCMS3/MSI/mucinous/inflammatory tumor program high
  poor normal-reference mapping

ambiguous:
  low UMI, high mitochondrial, doublet-like, only one weak tumor marker, or conflicting evidence
```

---

## 5. 肿瘤上皮主亚型：iCMS2 / iCMS3

### 5.1 iCMS2：canonical / WNT-MYC / proliferative / CIN-like

**核心程序：**

```text
WNT/TCF: ASCL2, LGR5, AXIN2, RNF43, NOTUM, NKD1, TCF7, LEF1
MYC/E2F/cell cycle: MYC, E2F1, E2F2, MKI67, TOP2A, MCM2-7, PCNA, TYMS, UBE2C
EGFR ligand axis: EREG, AREG
CIN support: chr8q/13q/20q gains, chr17p/18q losses
```

iCMS2 与 CMS2/canonical/WNT-MYC 更一致，Nature Genetics 2022 研究显示 iCMS2 上调 MYC/E2F targets，WNT 通路更强，并与典型染色体不稳定改变相关。([Nature][1])

**常见功能状态：**

```text
stem-like / crypt-base
transit-amplifying / proliferative
EGFR-ligand / canonical
CRIS-D-like WNT/IGF2
PDS1-like canonical stem-proliferative
部分 PDS3-like differentiated slow-cycling cells
```

### 5.2 iCMS3：secretory / mucinous / metabolic / inflammatory / EMT-prone

**核心程序：**

```text
secretory/goblet: MUC2, TFF3, SPINK4, CLCA1, FCGBP, AGR2, REG4
metabolic: SLC2A1, LDHA, HK2, PGK1, FABP1, APOA1, APOA4, MTTP, ALDOB
inflammatory/IFN: IFIT1, IFIT2, IFIT3, ISG15, MX1, OAS1, STAT1, IRF1, CXCL10, HLA-A/B/C, B2M
EMT/invasive support: LAMC2, LAMB3, ITGA5, MMP7, MMP14, SERPINE1, TGFBI, EMP1, L1CAM
```

iCMS3 在 bulk 上更常见于 CMS1/CMS3，常与 MSI-H、右侧、黏液性病理相关；iCMS3 上调 EMT、炎症、代谢相关通路，且 KRAS/PIK3CA 相对更常见。([Nature][1])

**常见功能状态：**

```text
goblet / mucinous / secretory
inflammatory / interferon-response
metabolic
EMT / invasive / tumor-budding
fetal / regenerative / oncofetal
hypoxia / stress
drug-tolerant / persister
```

---

## 6. 亚型内功能类别、marker 与判定逻辑

建议每个 malignant epithelial cell 同时给出：

```text
iCMS label + primary functional state + optional secondary state
```

例如：

```text
iCMS2; TA-proliferative
iCMS3; goblet-mucinous
iCMS3; EMT-invasive + hypoxia
iCMS2/iCMS3-mixed; fetal-regenerative
```

### 6.1 功能状态表

| 功能类别                                      | 主要 marker / signature                                                                                         | 常见亚型关联                                                             | 判定标准                                                               |
| ----------------------------------------- | ------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------ | ------------------------------------------------------------------ |
| **Stem-like / crypt-base**                | LGR5, OLFM4, ASCL2, SMOC2, AXIN2, RNF43, LRIG1, PROM1, CD44, SOX9                                             | iCMS2, CMS2-like, CRIS-D, PDS1                                     | WNT/TCF score 高；非正常 crypt stem；需 CNV/SNV 或 tumor-cluster 支持        |
| **Transit-amplifying / proliferative**    | MKI67, TOP2A, PCNA, MCM2-7, TYMS, CENPF, UBE2C, CDK1, CCNB1                                                   | iCMS2, PDS1                                                        | S/G2M score 高；不要让 cell cycle 单独驱动全部聚类                              |
| **Enterocyte-like / absorptive**          | CA1, CA2, AQP8, SLC26A3, CEACAM7, MS4A12, GUCA2A/B, KRT20, FABP1, APOA1/4                                     | iCMS3 或 PDS3-like differentiated；也可见 iCMS2 内 differentiated subset | 与正常 colonocyte 极易混淆；必须有 malignant 证据                               |
| **Goblet / mucinous / secretory**         | MUC2, TFF3, SPINK4, CLCA1, FCGBP, AGR2, REG4, MUC5AC, MUC5B                                                   | iCMS3, CMS3-like, CRIS-A                                           | 多个 mucin/secretory genes 同时高；排除 ambient mucin RNA                  |
| **BEST4+ absorptive-like**                | BEST4, OTOP2, SPIB, GUCA2B, CA7, NPY                                                                          | 多数为正常 mature absorptive；少数肿瘤 mimic                                 | 默认 normal-like；只有 CNV/SNV/tumor cluster 强阳性才标 malignant BEST4-like |
| **Enteroendocrine-like**                  | CHGA, CHGB, NEUROD1, PAX6, PYY, GCG, TPH1, SCG5                                                               | 少见；可为 lineage mimic 或 neuroendocrine differentiation               | 若 SYP/INSM1/DLL3/ASCL1 高，需单独考虑 neuroendocrine 成分                   |
| **Tuft-like**                             | DCLK1, TRPM5, POU2F3, AVIL, SH2D6, GFI1B                                                                      | 少见；正常 tuft 或 tumor tuft-like                                       | 需排除正常 tuft；结合 CNV/SNV                                              |
| **Paneth-like / deep crypt secretory**    | DEFA5, DEFA6, LYZ, REG3A, PLA2G2A, MMP7, KIT, SOX9                                                            | CRIS-E-like；右侧/IBD/化生背景更常见                                         | 不能只凭 LYZ；结合解剖部位和病理                                                 |
| **Inflammatory / IFN-response**           | IFIT1/2/3, ISG15, MX1, OAS1/2, STAT1, IRF1, CXCL9/10/11, HLA-A/B/C, B2M, TAP1                                 | iCMS3, CMS1-like, MSI-like                                         | 肿瘤上皮自身 ISG/MHC 高；排除 PTPRC+ 免疫 doublet                              |
| **Metabolic**                             | glycolysis: SLC2A1, LDHA, HK2, PGK1；absorptive/lipid: FABP1, APOA1/4, MTTP, ALDOB；mitochondrial/OXPHOS module | iCMS3, CMS3-like, CRIS-A                                           | 用 pathway score，不靠单基因；区分 hypoxia-driven glycolysis                 |
| **EMT / invasive / tumor-budding / pEMT** | LAMC2, LAMB3, ITGA5, ITGB1, MMP7, MMP14, PLAU, SERPINE1, TGFBI, EMP1, L1CAM, VIM, S100A4, ZEB1/2, SNAI2       | iCMS3, CRIS-B, CMS4/F overlay                                      | 必须仍有 epithelial/malignant 证据；排除 CAF                                |
| **Fetal / regenerative / oncofetal**      | CLU, ANXA1, L1CAM, TACSTD2, KRT17, KRT19, SOX9, EMP1, AREG, EREG, FOS, JUN, YAP/TEAD targets                  | PDS2-like, iCMS3-like, therapy/plasticity states                   | 常与损伤修复、复发、治疗耐受相关；可加 spatial/therapy 信息                             |
| **Hypoxia / stress / UPR**                | CA9, VEGFA, SLC2A1, LDHA, ENO1, PGK1, NDRG1, BNIP3, ADM, HSPA1A/B, DNAJB1, ATF3, DDIT3, XBP1, HSPA5           | 任意 iCMS；侵袭边缘/坏死区常见                                                 | 若所有细胞类型同时高，考虑 dissociation/stress artifact                         |
| **Drug-tolerant / persister**             | MEX3A, EMP1, L1CAM, CLU, ANXA1, KRT17, ALDH genes, ABCB1/ABCG2；低 MKI67/低 LGR5-WNT                             | post-treatment, MAPK/EGFR/chemo pressure                           | 需结合治疗背景；单样本未治疗数据中只能标 “DTP-like”                                    |

MEX3A+ 细胞被报道为 CRC drug-tolerant persister 状态，具有低 canonical LGR5/WNT 程序和再生样特征；近年研究也强调 oncofetal/YAP/AP-1 或 MAPK 驱动的 epithelial plasticity 与治疗耐受、再生样状态相关。([PubMed][6])

---

## 7. 经典 CMS、iCMS、CRIS、IMF 的联合判定

### 7.1 CMS 只做 sample / region pseudobulk

对每个样本或空间区域，把 malignant epithelial、immune、stromal 分开计算：

```text
CMS1-like:
  MSI/MMR-deficient, high immune infiltration, IFN/MHC, T cell/B cell/myeloid high
  epithelial correlate often iCMS3 + inflammatory/IFN

CMS2-like:
  WNT/MYC/canonical, CIN, epithelial-rich
  epithelial correlate often iCMS2 + stem/TA/proliferative

CMS3-like:
  metabolic, KRAS, mucinous/secretory tendency
  epithelial correlate often iCMS3 + goblet/metabolic

CMS4-like:
  stromal-rich, TGFβ, CAF/ECM, angiogenesis, invasion
  epithelial correlate can be pEMT/invasive, but CMS4 is mainly TME/fibrosis overlay
```

**关键规则：**

```text
不要把一个 CAF-rich 样本里的所有上皮细胞都叫 CMS4。
CMS4-like 应写成：sample-level CMS4/fibrotic; malignant epithelial cells include EMT-invasive/pEMT state.
```

### 7.2 IMF 标签

IMF 是对 iCMS 体系的修正，推荐在 sample/region 级别输出：

```text
iCMS2_MSS_NF
iCMS2_MSS_F
iCMS3_MSS_NF
iCMS3_MSS_F
iCMS3_MSI
```

判定逻辑：

```text
I = epithelial intrinsic state:
  iCMS2 or iCMS3 from malignant epithelial pseudobulk

M = MSI:
  clinical MSI/MMR IHC/WES/TMB > scRNA expression
  MLH1/MSH2/MSH6/PMS2 loss, MSI-H, high mutation burden

F = fibrosis:
  CAF/ECM/endothelial/spatial pathology support
  COL1A1/COL1A2/COL3A1/DCN/LUM/FAP/POSTN/THY1/ACTA2/TAGLN high
  TGFβ/ECM remodeling high
```

### 7.3 CRIS 辅助标签

```text
CRIS-A: mucinous / glycolytic / MSI or KRAS tendency
CRIS-B: TGFβ / EMT / invasive / poor prognosis
CRIS-C: EGFR signaling / EGFR inhibitor sensitivity tendency
CRIS-D: WNT activation / IGF2
CRIS-E: Paneth-like / TP53 tendency
```

实际 scRNA 中建议：

```text
per-patient malignant epithelial pseudobulk → CRIS classifier
clone/cluster pseudobulk → CRIS score
single-cell → 只报告 CRIS-like score，不硬分
```

---

## 8. 最新 atlas / 预印本如何纳入

大队列 CRC atlas 已经从“手工 marker 注释”转向 reference mapping。ICBI CRC atlas 整合超过 4.27 million cells、650 patients、49 studies，并提供 h5ad 和 scArches model 用于新数据映射；Broad/HTAN 相关 colon cancer atlas 也整合了 MMRp/MMRd 肿瘤、邻近正常和空间数据，用于识别 epithelial/TME programs。([crc.icbi.at][7])

2026 bioRxiv “perturbation-guided mapping” 方向提出 comparative single-cell CRC atlas 与 perturbation atlas 联合，强调超过 300 patients、1.5 million cells 的 epithelial atlas、hybrid/endoderm-like/oncofetal malignant states，以及用扰动数据把 cell state 映射到潜在机制；其中 endoderm-like/hybrid states 被描述为在 MSS KRAS-mutant CRC 中富集，并出现 ASCL2、EREG/AREG、HNF1A/FOXA3/GATA6/HNF4A/TCF4/SOX4/PROX1 等程序。该方向应作为扩展层使用，因为目前属于预印本，尚未同行评议。([ResearchGate][8])

落地方式：

```text
1. 先用本地规则完成 epithelial / malignant / iCMS / functional state 注释。
2. 再把 malignant epithelial cells 投影到 CRC atlas reference。
3. 对照 reference label：canonical, goblet, absorptive, inflammatory, EMT, fetal/regenerative, hybrid-endoderm-like 等。
4. 若 reference label 与本地 marker 冲突，以 CNV/SNV + marker + spatial/pathology 证据综合判定。
5. 预印本 label 只作为 secondary annotation，例如:
   "fetal-regenerative; Epi-CRC hybrid-endoderm-like score high; provisional"
```

空间转录组和多重成像可用于校正 scRNA 的 dissociation bias：例如 invasive margin、tumor buds、mucin pools、免疫排斥区域、fibrotic niche 等不能只靠单细胞表达推断；空间 atlas 研究显示 CRC 中 epithelial programs、immune exclusion、morphological gradients 与组织结构密切相关。([Harvard Tissue Atlas][9])

---

## 9. 推荐阈值与判定逻辑

### 9.1 Signature scoring

推荐使用：

```text
UCell / AUCell / singscore / AddModuleScore / ssGSEA
```

对每个样本内的 malignant epithelial cells 标准化：

```text
z_score = (cell_score - sample_mean) / sample_sd
```

### 9.2 单细胞功能状态阈值

建议规则：

```text
primary_state:
  最高功能 signature score ≥ 75th percentile
  且比第二高 signature 高 ≥ 0.15 z-score
  且 marker genes 至少 3–5 个共同表达

secondary_state:
  第二高 signature 也 ≥ 70th percentile
  且生物学兼容，例如 EMT+hypoxia, IFN+MSI, goblet+metabolic

mixed/ambiguous:
  top1-top2 < 0.15 z-score
  或只有单个 marker 高
  或该状态在所有细胞类型中普遍升高，提示 stress/ambient artifact
```

### 9.3 iCMS2/iCMS3 阈值

```text
iCMS2:
  iCMS2_score - iCMS3_score ≥ 0.15–0.25 z-score
  WNT/MYC/E2F module 高
  常伴 CIN/CNV evidence

iCMS3:
  iCMS3_score - iCMS2_score ≥ 0.15–0.25 z-score
  secretory/metabolic/inflammatory/EMT module 高
  可有 MSI/右侧/黏液性/KRAS/PIK3CA 支持

iCMS-mixed:
  两者差异小
  或同一 clone 内明显分化轨迹连续

iCMS-indeterminate:
  UMI 低、stress 高、doublet 风险高、marker 不足
```

对 patient/sample 层面：

```text
若 >60% malignant epithelial cells 为同一 iCMS → sample iCMS2 或 iCMS3
若 40–60% 混合 → report intratumoral heterogeneity
若空间区域差异明显 → region-level iCMS/IMF
```

### 9.4 PDS/SMI

PDS classifier 文献使用 0.6 默认概率阈值，低于阈值可标 mixed。实际 scRNA 中建议：

```text
PDS1-like:
  LGR5/ASCL2/MYC/E2F/cell-cycle 高

PDS2-like:
  ANXA1/REG regenerative + stromal/immune context 高

PDS3-like:
  slow-cycling, differentiated absorptive/EEC-like, PRC-high/MYC-low
```

PDS/SMI 尤其适合解释 **iCMS2 内部差异**，不要把它当作 iCMS 的替代品。([Nature][3])

---

## 10. 最常见混淆点与排除标准

| 混淆点                                  | 解决方案                                                                                |
| ------------------------------------ | ----------------------------------------------------------------------------------- |
| 正常 stem/TA 被误标为肿瘤 stem/proliferative | 必须有 CNV/SNV 或 patient-specific tumor cluster；正常 crypt stem/TA 在邻近正常中也会出现            |
| iCMS3/MSI 细胞 CNV 弱，被误判正常             | 用 MSI/MMR、突变、患者特异聚类、肿瘤程序和正常 reference mapping 联合判定                                  |
| EMT 细胞与 CAF 混淆                       | EMT 肿瘤细胞应保留 KRT8/18/19 或 EPCAM/CDH1 残留，并有 CNV/SNV；CAF 高 COL1A1/DCN/LUM/PDGFRA/ACTA2 |
| CMS4 被误当作单细胞上皮亚型                     | CMS4 多为 stromal/fibrotic bulk 信号；上皮细胞只能标 EMT-invasive/pEMT，样本标 CMS4-like/F          |
| IFN/inflammatory 被免疫 doublet 污染      | 排除 PTPRC/CD3D/LST1/MS4A1 高细胞；检查 HLA/ISG 是否由上皮自身表达                                   |
| Goblet/mucinous 被 ambient RNA 影响     | 要求 MUC2/TFF3/SPINK4/CLCA1/FCGBP 多 marker 一致，而不是单个 MUC2 高                            |
| BEST4+、absorptive、goblet 正常谱系被误标恶性   | 这些谱系默认先按 normal-like 处理，只有 CNV/SNV/tumor cluster 支持才标 malignant                     |
| Paneth-like 过度解释                     | 结肠中 Paneth-like/deep crypt secretory 需结合右侧、IBD、化生、肿瘤背景；LYZ 单独不够                     |
| Hypoxia/stress 是消化伪影                 | 若 HSPA/FOS/JUN/MT genes 在所有细胞类型同时升高，应标 dissociation/stress artifact                 |
| per-cell CMS/CRIS 不稳定                | CMS/CRIS 用 pseudobulk/cluster/sample；单细胞只报告 score 或 “like”                          |

---

## 11. 一套可执行的注释流程

```text
Step 0. QC
  remove low UMI, high mitochondrial, high ribosomal artifacts, doublets
  run SoupX/DecontX if mucin/ambient RNA severe

Step 1. Major lineage
  epithelial = EPCAM/KRT8/KRT18/KRT19/CDH1/CLDN4 high
  remove immune/stromal/endothelial/pericyte contaminants

Step 2. Normal vs malignant epithelial
  run CNV inference with patient-matched non-epithelial reference
  overlay patient-specific mutations
  map to matched normal + public gut/CRC atlas
  assign malignant-high / malignant-probable / normal-epithelial / ambiguous

Step 3. Tumor epithelial intrinsic subtype
  score iCMS2 and iCMS3
  assign iCMS2/iCMS3/mixed/indeterminate
  perform patient/clone pseudobulk confirmation

Step 4. Functional state annotation
  score stem, TA, absorptive, goblet, BEST4, EEC, tuft, Paneth-like,
        IFN, metabolic, EMT-invasive, fetal-regenerative, hypoxia, DTP
  assign primary and secondary functional states

Step 5. Orthogonal classifiers
  sample/region pseudobulk CMS
  cluster/pseudobulk CRIS-A–E
  patient/region IMF = iCMS + MSI + fibrosis
  PDS/SMI for CMS2/iCMS2 internal structure

Step 6. Spatial / multi-omics correction
  validate EMT/invasive at tumor front
  validate hypoxia near necrosis
  validate F/fibrosis with CAF/ECM-rich regions
  validate mucinous/goblet with pathology
  validate DTP/regenerative with treatment history

Step 7. Final report
  per-cell labels + confidence
  per-cluster marker table
  per-patient subtype proportions
  ambiguous/excluded populations listed separately
```

---

## 12. 建议最终命名格式

推荐命名：

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

最重要的执行准则是：**先证明“恶性上皮”，再标 iCMS2/iCMS3，最后标功能状态**。CMS、CRIS、IMF、PDS、Epi-CRC atlas 都应作为正交解释层，而不是彼此互相替代。
