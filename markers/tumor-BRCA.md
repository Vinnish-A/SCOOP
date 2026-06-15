
**第一步：先判定大类细胞。**
上皮细胞用 `EPCAM, TACSTD2, KRT8, KRT18, KRT19, CDH1, MUC1`；排除免疫 `PTPRC`、内皮 `PECAM1, VWF, KDR`、成纤维/基质 `COL1A1, DCN, LUM`、周细胞/平滑肌 `RGS5, MCAM, ACTA2`、脂肪细胞 `ADIPOQ, PLIN1`。

**第二步：在上皮细胞内区分肿瘤与正常。**
肿瘤上皮细胞不能只靠 EPCAM 或 KRT 表达判断，应优先用：

1. scRNA 推断 CNV：inferCNV、CopyKAT、CONICSmat 或 scCNA；
2. 与患者肿瘤 WES/scDNA/FISH/IHC 的突变或扩增一致性；
3. 与患者内克隆结构一致的 arm-level CNA；
4. 显著偏离 Human Breast Cell Atlas 正常上皮参考；
5. 空间上位于肿瘤巢而非正常导管/小叶双层结构。

Pang 等在乳腺癌 scRNA 整合分析中使用 inferCNV 识别 malignant cells，并强调 CNV 亚克隆与转录状态并非一一对应，因此 **CNV 用于判定 malignant，不应直接等同于功能状态**。([Nature][2])

**第三步：正常上皮先单独注释，再只对 malignant epithelial cells 做肿瘤亚型和功能状态。**
Human Breast Cell Atlas 已经提供健康乳腺参考：Reed 等报道 55 名供体、80 万以上细胞；Kumar 等构建了空间分辨的成人乳腺图谱；Bhat-Nakshatri 等 snRNA/snATAC 图谱给出了 LHS、LASP、basal-myoepithelial 的核心身份基因；iHBCA v1.0 进一步整合 287 名供体、约 212 万细胞。([Nature][3])

## 3. 正常乳腺实质上皮细胞注释

| 正常上皮类别                                          | 推荐标签                                          | 核心 marker                                                         | 判断标准                                                               |
| ----------------------------------------------- | --------------------------------------------- | ----------------------------------------------------------------- | ------------------------------------------------------------------ |
| Luminal epithelial, general                     | `NormalEpi_Luminal`                           | `EPCAM, KRT8, KRT18, KRT19, CDH1, MUC1`                           | CNV-neutral；与正常乳腺 luminal reference 接近；无患者特异性肿瘤 CNA/突变             |
| Hormone-sensing luminal cell                    | `NormalEpi_LHS`                               | `ESR1, PGR, AR, FOXA1, GATA3, TFF1, TFF3, BCL2`                   | luminal core 阳性；`ESR1/PGR/FOXA1` 高；`MKI67/TOP2A` 低；CNV-neutral     |
| Luminal progenitor / LASP / secretory precursor | `NormalEpi_LASP`                              | `EHF, ELF5, KIT, ALDH1A3, LTF, SLPI, S100A14, KRT8, KRT18, KRT19` | luminal core 阳性；`ESR1/PGR` 低或阴性；`EHF/ELF5` 高；CNV-neutral           |
| Mature secretory / alveolar / lactocyte-like    | `NormalEpi_Secretory` 或 `NormalEpi_Lactocyte` | `LALBA, CSN2, CSN3, LTF, WFDC2, ELF5, MUC1`                       | 妊娠/哺乳或局部 secretory 程序背景下更可信；需 CNV-neutral                          |
| Basal / myoepithelial cell                      | `NormalEpi_BasalMyo`                          | `KRT5, KRT14, KRT15, KRT17, TP63, ACTA2, MYL9, TAGLN, ITGA6`      | basal/myoepithelial marker 高；常与 luminal layer 构成正常双层结构；CNV-neutral |
| Cycling epithelial cell                         | `NormalEpi_Cycling_{parent}`                  | `MKI67, TOP2A, UBE2C, CDC20, CDK1, CCNB1, PCNA, MCM2/5/6`         | 必须同时有正常 parent marker 且 CNV-neutral；不要因为 Ki-67 高就自动判为肿瘤            |

实践中，`Cycling` 应作为附加状态，而不是独立细胞类型。例如 `NormalEpi_Cycling_LASP` 比 `Cycling epithelial` 更精确。2025 年乳腺上皮命名共识也强调，目前高置信正常乳腺上皮命名应围绕 basal-myoepithelial、luminal hormone-sensing、luminal adaptive secretory precursor 等清晰类别。([Cell][4])

## 4. 肿瘤细胞经典亚型注释

只在 `MaligEpi` 中做以下亚型判断。建议使用 `PAM50/SCSubtype score + marker + CNV/IHC prior`，而不是单基因判定。Wu 等的乳腺癌单细胞与空间图谱提出了 SCSubtype，用于在单细胞层面识别 recurrent neoplastic cell intrinsic subtype heterogeneity；Xu 等 2024 年整合 236,363 个细胞、119 个乳腺肿瘤样本，并在分子亚型和 10 类上皮转录异质性层面解析癌上皮细胞。([单细胞门户][5])

| 肿瘤亚型                           | 推荐标签                                                 | 正向 marker / signature                                                                      | 排除或低表达                                                   | 判断标准                                                                                                                                |
| ------------------------------ | ---------------------------------------------------- | ------------------------------------------------------------------------------------------ | -------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| Luminal A-like                 | `MaligEpi_LuminalA`                                  | `ESR1, PGR, FOXA1, GATA3, XBP1, TFF1, TFF3, BCL2, PIP, SCGB2A2, AGR2/3`                    | `MKI67, TOP2A, UBE2C, CCNB1` 低；`ERBB2` 低；basal markers 低 | CNV+；HR/luminal score 高；proliferation score 低；临床 ER/PR+、HER2- 更支持                                                                   |
| Luminal B-like                 | `MaligEpi_LuminalB`                                  | Luminal markers 同上，同时 `MKI67, TOP2A, UBE2C, CDC20, CDK1, CCNB1/2, AURKB` 高                 | basal markers 低；HER2 可低/中/高                              | CNV+；HR/luminal score 高；proliferation score 高；常表现为 ER+ 但 Ki-67 高或 PR 相对低                                                            |
| HER2-enriched                  | `MaligEpi_HER2E`                                     | `ERBB2, GRB7, STARD3, MIEN1, PGAP3, TCAP, KRT7, KRT19`；可见 17q12 amp                        | `ESR1/PGR` 常低或中等；不要求绝对阴性                                 | CNV+；HER2/ERBB2 module 高；ERBB2 amp 或 IHC/ISH HER2+ 强支持；若 ER 高同时 HER2 高，标为 `Luminal_HER2activated` 更稳妥                               |
| Basal-like / TNBC-like         | `MaligEpi_BasalLike` 或 `MaligEpi_TNBC_BasalLike`     | `KRT5, KRT14, KRT17, KRT6A/B/C, EGFR, TP63, ITGA6, LAMC2, S100A2, CD44`，常伴 `MKI67/TOP2A` 高 | `ESR1, PGR, ERBB2` 低                                     | CNV+；basal score 高；HR/HER2 低。临床 TNBC 是 ER/PR/HER2 阴性定义，Basal-like 是转录亚型，二者高度相关但不完全等价                                                |
| Claudin-low / mesenchymal-like | `MaligEpi_ClaudinLowLike` 或作为 tag：`EMThi-Claudinlow` | `VIM, FN1, ZEB1/2, SNAI2, TWIST1, AXL, CDH2, COL1A1/2, ACTA2, SPARC, ITGA5, CD44, ALDH1A1` | `CLDN3, CLDN4, CLDN7, EPCAM, CDH1, MUC1, CD24` 低         | 建议作为跨亚型 phenotype/tag，而非强行互斥主亚型；Prat 等描述 claudin-low 为低 tight-junction、EMT/stem/immune 富集，后续研究也提示它更像可横跨亚型的 mesenchymal/stem-like 状态 |
| TNBC-LAR / AR-luminal-like     | `MaligEpi_TNBC_LAR` 或 `BasalLike_ARhi`               | `AR, FOXA1, GATA3, XBP1, KRT8, KRT18, MUC1`；ER 通常低                                         | `ESR1/PGR/ERBB2` 低；basal markers 可不高                     | 只在临床 TNBC 或 HR/HER2 低的 malignant cells 中使用；TNBCtype-4 将 TNBC 分为 BL1、BL2、M、LAR，LAR 为肿瘤内在 AR/luminal-like 亚型                          |
| Normal-like                    | 不建议作为肿瘤主标签                                           | 正常 luminal/basal/stromal admixture signature                                               | —                                                        | 单细胞中优先解释为正常上皮/低肿瘤纯度/双细胞或 ambient RNA；只有 CNV+ 且确有肿瘤克隆证据时才可保留 `MaligEpi_NormalLike-like`                                              |

Claudin-low 的原始定义包括低 luminal marker、低 claudin/黏附分子、高 EMT/immune/stem-like 特征，多数为 ER-/PR-/HER2-；但近年更推荐把它看成一种可横跨 intrinsic subtype 的 phenotype，而不是绝对独立的 PAM50 同级类别。([PubMed][6]) TNBC 方面，Lehmann 的 TNBCtype-4 使用 BL1、BL2、M、LAR 四类肿瘤内在亚型；多组学研究也确认 TNBC 具有 basal-like、mesenchymal、LAR 等不同生物学状态。([PLOS][7])

## 5. 亚型内功能状态：推荐统一功能标签

这些功能标签**可叠加**。例如一个细胞可以是 `LuminalB_HRhi-Prolifhi-Stresshi`，也可以是 `BasalLike_Prolifhi-IFN_MHCIIhi`。Pang 等 2024 年整合乳腺癌 scRNA-seq，提出跨患者复现的 consensus cancer cell states：hormone response、protein folding、lipid metabolism/HER2-like、proliferation、immune、hypoxia、EMT 等；其 hc3 高表达 `CDC20, CDK1, MKI67`，hc10 高表达 `VIM, COL1A1, ACTA2`，hc7 具有 MHC-II/CD74 相关免疫交互特征。([Nature][2])

| 功能状态                                            | 推荐 tag                           | marker                                                                                                               | 适用亚型                                 | 判定标准                                                                                                                                              |
| ----------------------------------------------- | -------------------------------- | -------------------------------------------------------------------------------------------------------------------- | ------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| Hormone response / ER signaling                 | `HRhi`                           | `ESR1, PGR, FOXA1, GATA3, XBP1, TFF1, TFF3, BCL2, PIP, SCGB2A2`                                                      | Luminal A/B，部分 HER2+ ER+             | HR module 位于 malignant cells 前 25–30%；至少 3 个核心基因表达；`ESR1/PGR` dropout 时看 FOXA1/GATA3/XBP1/TFF1                                                    |
| Proliferative / Cycling                         | `Prolifhi`                       | `MKI67, TOP2A, UBE2C, CDC20, CDK1, CCNB1/2, AURKB, MCM2/5/6, PCNA`                                                   | Luminal B、HER2E、Basal-like 常见        | G2M/S score 高；cluster 中 ≥30–40% 细胞高表达；用于区分 Luminal A vs Luminal B                                                                                 |
| HER2 / RTK activated                            | `HER2Signalhi`                   | `ERBB2, GRB7, STARD3, PGAP3, MIEN1, EGFR, ERBB3, AREG`                                                               | HER2E、Luminal-HER2 activated、部分 TNBC | ERBB2/GRB7 module 高；若有 17q12 amp/IHC 3+ 或 ISH+ 证据则高置信                                                                                             |
| Basal-stem / squamous-like                      | `BasalStemhi`                    | `KRT5, KRT14, KRT17, KRT6A/B/C, TP63, EGFR, ITGA6, LAMC2, S100A2, CD44`                                              | Basal-like/TNBC，部分 metaplastic-like  | basal score 高且 HR/HER2 低；若同时 `EPCAM/CDH1` 保留，可标 partial basal                                                                                     |
| Claudin-low / EMT / mesenchymal                 | `EMThi` 或 `ClaudinLowhi`         | 高：`VIM, FN1, ZEB1/2, SNAI2, TWIST1, AXL, CDH2, COL1A1/2, ACTA2, SPARC, SDC2`；低：`CLDN3/4/7, EPCAM, CDH1`              | Claudin-low-like、TNBC-M、侵袭前沿         | EMT module 高；junction/epithelial module 低；需排除 CAF doublet，因为 `COL1A1/ACTA2` 也来自基质                                                                 |
| Hybrid epithelial–mesenchymal plasticity        | `HybridEMP`                      | 同时表达 `EPCAM/CDH1/KRT8/18` 与 `VIM/FN1/S100A2/CRYAB/KRT14`                                                             | 转移、侵袭边缘、治疗后样本                        | epithelial 与 mesenchymal module 均高；不要仅因 VIM 单基因阳性判定。乳腺癌转移单细胞研究显示 EMP 是 ITH 的主要来源，并有中间 EMP 状态与较差结局相关。([JCI Insight][8])                            |
| Hypoxia / glycolysis                            | `Hypoxiahi`                      | `CA9, VEGFA, HILPDA, NDRG1, BNIP3, LDHA, SLC2A1, ENO1, PGK1`                                                         | 各亚型均可出现，TNBC/HER2E 常见                | hypoxia/HIF1A target score 高；空间上靠坏死/低氧区域更可信                                                                                                       |
| Stress / UPR / protein folding                  | `StressUPRhi`                    | `HSPA1A/B, HSP90AA1, HSPA5, DDIT3, ATF4, XBP1, DNAJB1, HSPB1`                                                        | 治疗后、低质量/高应激样本、ER+ hc1-like           | 同时检查 QC：若线粒体比例高、UMI 低，标为 `Stress_artifact_possible`                                                                                               |
| Antigen presentation / IFN / immune-interacting | `IFN_MHChi` 或 `ImmuneInteracthi` | `HLA-A/B/C, B2M, TAP1/2, NLRC5, HLA-DRA/DRB1/DPA1/DPB1, CD74, STAT1, IRF1, ISG15, IFI6, IFI27, IFIT1/3, CXCL9/10/11` | TNBC、HER2E、部分 ER+ 炎症样区域              | 必须同时有 epithelial marker 与 CNV+，否则易误判为免疫 doublet；Xu 等 InteractPrint 也强调癌上皮异质性会影响 T cell interaction 和 anti-PD-1 response。([OA Monitor Ireland][9]) |
| Secretory / lactation-like                      | `Secretoryhi`                    | `PIP, SCGB2A2, LTF, SLPI, WFDC2, ELF5, AGR2/3, MUC1`                                                                 | Luminal tumor 或正常 LASP/secretory     | 关键是 CNV：CNV-neutral 多为正常 LASP/secretory；CNV+ 才作为肿瘤分泌样状态                                                                                           |
| AR / luminal androgen receptor                  | `ARhi` 或 `LAR`                   | `AR, FOXA1, GATA3, XBP1, KRT8, KRT18, MUC1`                                                                          | TNBC-LAR 或 ER-low luminal-like       | 在 `ESR1/PGR/ERBB2` 低背景中 AR module 高；不要把普通 ER+ luminal 误标为 LAR                                                                                     |
| Inflammatory cytokine / chemokine               | `InflamChemohi`                  | `CXCL1, CXCL2, CXCL3, CXCL8, IL6, CCL2, CSF1, CXCL16`                                                                | TNBC、治疗后、坏死/炎症区域                     | 与 IFN/MHC 分开；需排除髓系 doublet 和 ambient RNA                                                                                                          |
