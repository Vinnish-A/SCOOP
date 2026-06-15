
## 2. 肿瘤细胞经典亚型与亚型内功能类别

### 2.1 AC-like / astrocyte-like，常与 Classical/astrocytic 程序相关

**核心判定：** 恶性 CNA 阳性或突变阳性，且 AC-like module score 最高。

| 类别                                     | Marker / signature                                                  | 判断标准                                                                                                    |
| -------------------------------------- | ------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| AC-like 主状态                            | `GFAP`, `AQP4`, `ALDOC`, `APOE`, `SLC1A3`, `S100B`, `MLC1`, `SOX9`  | AC score 为四状态最高；同时有肿瘤 CNA/SNV；不能仅凭 GFAP/AQP4 判恶性，因为正常星形胶质细胞也表达                                          |
| AC-homeostatic / differentiated        | `AQP4`, `ALDOC`, `SLC1A2`, `SLC1A3`, `GLUL`, `GFAP`，低 `MKI67/TOP2A` | AC 主状态 + 低 cycling、低 hypoxia、低 stress                                                                   |
| AC-connectivity / tumor-microtube-like | `CHI3L1`, `GAP43`, `APOE`, `GJA1`, `TTYH1`，可加入空间或形态证据               | AC-like 或 MES-like 中连接性 signature 高；研究显示 AC-like/MES-like 连接性更高，`CHI3L1` 是网络连接和不良预后的重要标志物。([Nature][3]) |
| AC-gliosis / injury-response           | `VIM`, `SERPINA3`, `LCN2`, `C3`, `CD44`, `CHI3L1`, `GFAP`           | AC score 高 + gliosis/injury score 高；GBM-Space 预印本将传统 MES 重新解释为胶质损伤反应与缺氧响应连续谱的一部分。([ResearchGate][4])    |
| AC-cilia-like                          | `DNAH*`, `CFAP*`, `FOXJ1`, `TPPP3`                                  | CNA 阳性 + cilia module 高；CARE 2025 发现 cilia-like 是少量、约 1.6% 的 AC-like 相关恶性状态。([Nature][5])               |

**注意：** AC-like 与正常 astrocyte marker 高度重叠，所以 AC-like 标签必须以 CNA/SNV 或肿瘤参考映射为前提。

---

### 2.2 OPC-like / oligodendrocyte precursor-like，常与 Proneural/PDGFRA 程序相关

| 类别                         | Marker / signature                                                               | 判断标准                                                                                |
| -------------------------- | -------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| OPC-like 主状态               | `PDGFRA`, `OLIG1`, `OLIG2`, `SOX10`, `CSPG4`, `PTPRZ1`, `BCAN`, `VCAN`, `NKX2-2` | OPC score 最高 + CNA/SNV 阳性；PDGFRA 扩增支持 OPC-like                                      |
| OPC-progenitor / stem-like | `PDGFRA`, `CSPG4`, `PTPRZ1`, `OLIG2`, `SOX10`, `NES`, `PROM1`                    | OPC 主状态 + progenitor/stemness score 高                                               |
| OPC-cycling                | `MKI67`, `TOP2A`, `PCNA`, `MCM2-7`, `CENPF`, `UBE2C`, `CDK1`, `CCNB1/2`          | OPC 主状态 + S/G2M score 高；不要把 cycling 单独作为第五经典亚型，而作为功能 overlay                        |
| OPC-invasive / motile      | `TNC`, `ANXA1`, `ITGA6`, `ITGB1`, `MMP2`, `MMP14`, `SERPINE1`                    | OPC 主状态 + invasion/motility score 高；CARE 分析提示 OPC/NPC/NEU 相关程序更偏运动和侵袭。([Nature][5]) |
| OPC-oligo-differentiating  | `PLP1`, `MBP`, `MAG`, `MOG`, `MOBP`, `OPALIN`                                    | 只有 CNA/SNV 阳性时才可标为肿瘤性 oligo-differentiating；否则优先判为正常 oligodendrocyte                |

Couturier 等在 GBM 单细胞研究中提出 GBM 存在保守的神经三谱系层级，中心偏 glial progenitor-like，且多数 cycling cell 富集于 progenitor 程序中；这支持把 OPC/GPC/progenitor 与 cycling 拆开注释。([Nature][6])

---

### 2.3 NPC-like / neural progenitor-like，常与 Proneural/CDK4 程序相关

| 类别                                      | Marker / signature                                                                   | 判断标准                                                                                                                              |
| --------------------------------------- | ------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------- |
| NPC-like 主状态                            | `SOX2`, `SOX4`, `SOX11`, `ASCL1`, `DLL3`, `HES6`, `DCX`, `STMN2`, `TUBB3`, `NEUROD1` | NPC score 最高 + CNA/SNV 阳性；CDK4 扩增支持 NPC-like                                                                                      |
| NPC1-progenitor / Notch-developmental   | `SOX2`, `ASCL1`, `DLL3`, `HES6`, `SOX4`, `SOX11`                                     | GBmap 六状态下 NPC1 score 高；偏神经发育/干性                                                                                                  |
| NPC2-neuronal-lineage / neuroblast-like | `DCX`, `STMN2`, `TUBB3`, `MAP2`, `NEUROD1`                                           | GBmap 六状态下 NPC2 score 高；必须与正常 neuron 区分                                                                                           |
| NPC-cycling                             | `MKI67`, `TOP2A`, `MCM2-7`, `UBE2C`, `BIRC5`                                         | NPC 主状态 + cycling score 高                                                                                                         |
| NPC/NEU-like updated state              | `NRG1`, `NRG3`, `NRXN3`, `CNTNAP2`, `SYN3`, `SYT1`                                   | CARE 2025 新增 NEU-like malignant program；若 CNA 阳性且 NEU-like score 高，ClassicState 可归 NPC-like，UpdatedState 写 NEU-like。([Nature][5]) |

**注意：** `STMN2`, `SYT1`, `SNAP25`, `RBFOX3` 等可来自正常 neuron 或环境 RNA。若无 CNA/SNV，不能把 neuronal marker 阳性细胞标为 NPC/NEU-like 肿瘤细胞。

---

### 2.4 MES-like / mesenchymal-like，传统 Mesenchymal，但建议拆成 gliosis、hypoxia、stress、ECM 等功能

| 类别                               | Marker / signature                                                                                  | 判断标准                                                                                           |
| -------------------------------- | --------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------- |
| MES-like 主状态                     | `CD44`, `CHI3L1`, `VIM`, `LGALS3`, `S100A10`, `ANXA1`, `SERPINE1`, `TNC`, `ITGA5`, `STAT3`, `CEBPB` | MES score 最高 + CNA/SNV 阳性；NF1 loss/mutation 支持 MES-like                                        |
| MES1-inflammatory / gliosis-like | `CD44`, `CHI3L1`, `LGALS3`, `S100A10`, `ANXA1`, `VIM`, `CCL2`, `CXCL8`, `IL6`, `STAT3`, `CEBPB`     | MES score 高 + inflammatory/gliosis score 高；通常与髓系细胞丰富生态位相关                                      |
| MES2-hypoxia / angiogenic        | `VEGFA`, `CA9`, `HILPDA`, `ADM`, `BNIP3`, `NDRG1`, `SLC2A1`, `LDHA`, `EGLN3`                        | hypoxia score 高；建议 UpdatedState 写 Hypoxia，而不是笼统 MES2；空间研究显示缺氧是 glioma 组织结构的重要组织者。([PubMed][7]) |
| MES-stress / UPR / heat-shock    | `HSPA1A`, `HSPA1B`, `HSP90AA1`, `DNAJB1`, `ATF3`, `DDIT3`, `XBP1`                                   | stress score 高；CARE 将 stress 从传统 MES-like 中拆出。([Nature][5])                                    |
| MES-ECM / invasive-remodeling    | `TNC`, `FN1`, `COL1A1`, `COL1A2`, `COL3A1`, `ITGA5`, `MMP2`, `MMP14`, `SERPINE1`                    | ECM/remodeling score 高 + CNA 阳性；若无 CNA，需警惕 fibroblast/perivascular cell contamination          |
| MES-connectivity                 | `CHI3L1`, `GAP43`, `APOE`                                                                           | MES 主状态 + connectivity signature 高；`CHI3L1` 与肿瘤微管连接网络功能相关。([Nature][3])                        |

---

## 3. 最新权威大团队图谱 / 预印本如何并入命名

### 3.1 GBmap：推荐作为参考映射层

GBmap 整合了 240 名患者、26 个数据集、超过 110 万细胞，建立 IDH-wildtype GBM 单细胞/空间参考图谱；其模型把恶性细胞映射到 AC-like、MES1-like、MES2-like、NPC1-like、NPC2-like、OPC-like 六个状态，并建议对低 posterior probability 的细胞标记低置信度。([PMC][8])

**落地做法：**

* Primary ClassicState：把 GBmap 六状态折叠为四状态

  * AC-like → AC-like
  * OPC-like → OPC-like
  * NPC1/NPC2 → NPC-like
  * MES1/MES2 → MES-like
* UpdatedState：保留六状态

  * `NPC1-like` 偏 progenitor/developmental
  * `NPC2-like` 偏 neuronal-lineage
  * `MES1-like` 偏 inflammatory/gliosis
  * `MES2-like` 偏 hypoxia/angiogenic

### 3.2 CARE / Nature Genetics 2025：推荐作为高分辨率 malignant-state 字典

CARE 2025 的大规模单核分析重新审视了 GBM 恶性状态：在原有 AC、OPC、NPC、MES 基础上，将 hypoxia、stress 从 MES-like 中拆分，并发现新的 `NEU-like`、`GPC-like`、`cilia-like` 程序；同时指出 cycling 更应作为叠加特征，而不是独立细胞身份。([Nature][5])

**推荐纳入 UpdatedState 的十类：**

1. `AC-like`
2. `OPC-like`
3. `NPC-like`
4. `MES-like`
5. `Hypoxia`
6. `Stress`
7. `NEU-like`
8. `GPC-like`
9. `Cilia-like`
10. `Cycling overlay`，注意 cycling 是 overlay，不是主身份

其中 `GPC-like` 的代表 marker 可用 `EGFR`, `ALK`, `MEIS1`, `MEOX2`, `ETV1`, `ELOVL2`；`NEU-like` 可用 `NRG1`, `NRG3`, `NRXN3`, `CNTNAP2`, `SYN3`, `SYT1`；`Cilia-like` 可用 `DNAH*`, `CFAP*`。CARE 还显示约 20% 恶性细胞呈 hybrid state，且部分杂合状态可能代表真实状态转变而非技术 doublet。([Nature][5])

### 3.3 GBM-Space / Wellcome LEAP 预印本：建议把“亚型”理解为空间轨迹上的区域状态

GBM-Space 由 Wellcome LEAP Delta Tissue 团队构建多模态、多区域单细胞和空间图谱，包含超过 100 万单核转录组、12 个肿瘤、多个空间采样位点；其预印本强调 GBM 异质性来自从发育样状态到 astrocyte-like、再到 gliosis/hypoxia response 的空间化轨迹。([GBM-Space][9])

**在注释上应体现为：**

* 不把 MES-like 简单理解为固定谱系，而标成：

  * `MES-like; gliosis/injury-response`
  * `MES-like; hypoxia-response`
  * `AC-like; AC-gliosis transition`
* 对空间数据，增加 `Niche` 字段：

  * hypoxic / necrotic niche
  * macrophage-rich inflammatory niche
  * proliferative progenitor niche
  * neuronal-interface / invasive margin niche
  * perivascular / ECM niche

---

## 4. 判断标准：推荐的可执行打分规则

### 4.1 恶性细胞识别

先做 tumor/normal，再做 tumor state。推荐证据优先级：

1. **CNA inference：最高优先级**
   用 inferCNV、CopyKAT、Numbat、HoneyBADGER、CaSpER 等从 sc/snRNA 推断大尺度 CNV。GBM 常见证据包括 chr7 gain、chr10 loss、EGFR amplification、CDKN2A/B deletion、PDGFRA amplification、CDK4/MDM2 amplification、PTEN 区域缺失等。Patel 等早期 GBM 单细胞研究已使用 scRNA 推断 CNV 区分肿瘤特异性群体；后续 GBM 单细胞分析也常结合 PTPRC、marker 和 CNA 定义肿瘤细胞。([Broad Institute Portals][10])
2. **SNV / fusion / EGFRvIII / allele-specific CNA：强证据**
   若有 WES/WGS/targeted panel 或 multiome，可用肿瘤突变直接确认。
3. **参考图谱映射：辅助证据**
   用 GBmap / CARE malignant-state reference 映射恶性细胞；正常脑细胞则用成人脑参考图谱，例如 Adult Human MTG / CellTypist，而不是只依赖 GBM 肿瘤图谱。([Research Square][11])
4. **marker 只能作为辅助**
   因为 AC-like 与 astrocyte、OPC-like 与 OPC、NPC/NEU-like 与 neuron marker 高度重叠。

### 4.2 经典亚型分配

对每个恶性细胞计算四个模块分数：AC、OPC、NPC、MES。可用 UCell、AUCell、AddModuleScore、GSVA、scVI/scANVI classifier。

**建议规则：**

* `ClassicState = argmax(AC, OPC, NPC, MES)`
* 高置信度：最高分超过背景阈值，且 top1 - top2 ≥ 0.10–0.15，或按 CARE 类似做法采用更宽的 hybrid 阈值；CARE 使用 top two state scores 差距来识别 hybrid 状态。([Nature][5])
* Hybrid：top1 与 top2 接近，且不是 doublet，例如 `AC/MES hybrid`、`OPC/NPC hybrid`
* Low-confidence：所有 state score 低，或 GBmap posterior probability < 0.5

### 4.3 功能模块分配

功能模块不应互斥。一个细胞可以是 `OPC-like + cycling + invasive`，也可以是 `MES-like + hypoxia + ECM-remodeling`。

推荐功能模块阈值：

* cluster-level：模块分数高于所有恶性细胞均值 + 1 SD，且 marker 在该 cluster ≥20–30% 细胞表达；
* cell-level：模块分数在恶性细胞 top 10–15%，或 AUCell/UCell 显著高；
* 若两个功能模块同时高，允许多标签；
* cycling、hypoxia、stress 必须作为 overlay，不覆盖 ClassicState。

---

## 5. 如何区分正常实质细胞与肿瘤细胞

GBM 单细胞里最容易误判的是：正常 astrocyte 被误标为 AC-like tumor，正常 OPC 被误标为 OPC-like tumor，正常 neuron 被误标为 NPC/NEU-like tumor，正常 oligodendrocyte 被误标为 oligo-differentiated tumor。核心原则是：**marker 相似时，以 CNA/SNV/肿瘤突变优先。**

### 5.1 正常实质细胞包括哪些

严格意义上的脑实质细胞主要包括 neuron、astrocyte、oligodendrocyte、OPC；靠近脑室或样本含室管膜时还可见 ependymal/ciliated cells。microglia 是脑驻留免疫细胞，常与实质一起分析，但建议单独归入 myeloid/microglia，而不是神经实质谱系。

| 正常细胞                 | Marker                                                                                                                           | 与肿瘤混淆点                         | 区分标准                                                                    |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------- | ------------------------------ | ----------------------------------------------------------------------- |
| Astrocyte            | `AQP4`, `ALDH1L1`, `SLC1A2`, `SLC1A3`, `GLUL`, `GJA1`, `GFAP`, `S100B`, `MLC1`                                                   | 与 AC-like tumor 重叠             | 正常 astrocyte 应 copy-neutral，无 EGFR/chr7/chr10 等肿瘤 CNA                   |
| OPC                  | `PDGFRA`, `CSPG4`, `VCAN`, `PTPRZ1`, `OLIG1`, `OLIG2`, `SOX10`                                                                   | 与 OPC-like tumor 重叠            | 正常 OPC copy-neutral；若 PDGFRA/chr7/chr10 等异常则偏肿瘤                         |
| Oligodendrocyte      | `MBP`, `PLP1`, `MOG`, `MAG`, `MOBP`, `OPALIN`, `CLDN11`                                                                          | 与 oligo-like tumor 或环境 RNA 混淆  | 高髓鞘 marker + copy-neutral → 正常；CNA+ 才考虑肿瘤性分化                            |
| Neuron               | pan-neuronal: `RBFOX3`, `SNAP25`, `SYT1`, `MAP2`, `TUBB3`; excitatory: `SLC17A7`, `CAMK2A`; inhibitory: `GAD1`, `GAD2`, `SLC6A1` | 与 NPC/NEU-like tumor、环境 RNA 混淆 | 正常 neuron 通常 copy-neutral、成熟神经 marker 强；CNA+ 且 NEU-like signature 高才标肿瘤 |
| Ependymal / ciliated | `FOXJ1`, `PIFO`, `TPPP3`, `DNAH5`, `CFAP*`                                                                                       | 与 CARE cilia-like tumor 混淆     | copy-neutral → 正常 ependymal/ciliated；CNA+ → malignant cilia-like        |
| Microglia            | `P2RY12`, `TMEM119`, `CX3CR1`, `SALL1`, `CSF1R`, `AIF1`                                                                          | 与 macrophage/TAM 混淆，不应归肿瘤      | `PTPRC` 阳性、copy-neutral；肿瘤细胞一般不应强表达 pan-immune marker                   |
| Endothelial          | `PECAM1`, `VWF`, `CLDN5`, `KDR`, `FLT1`                                                                                          | 血管生态位样本混入                      | copy-neutral + endothelial marker                                       |
| Pericyte / mural     | `PDGFRB`, `RGS5`, `MCAM`, `ACTA2`, `TAGLN`                                                                                       | 与 MES/ECM tumor 混淆             | copy-neutral；强 mural marker                                             |
| Fibroblast / VLMC    | `COL1A1`, `COL1A2`, `DCN`, `LUM`, `COL3A1`                                                                                       | 与 MES-ECM tumor 混淆             | 无 CNA → 非肿瘤基质；CNA+ 且 MES score 高 → MES-ECM tumor                        |

### 5.2 推荐判定流程

1. **先做 broad cell type annotation**
   用 normal brain reference 把 copy-neutral 细胞标为 neuron、astrocyte、OPC、oligodendrocyte、ependymal、microglia、endothelial、pericyte、fibroblast/VLMC、immune cells。

2. **再对疑似肿瘤群做 CNA/SNV 判定**

   * CNA+ / tumor SNV+：进入 GBM malignant annotation
   * CNA− / tumor SNV−：保留为正常或 TME
   * CNA 弱但 tumor-state score 高：标为 `putative malignant, low-confidence`，需空间或突变验证

3. **排除 doublet 与 ambient RNA**

   * `PTPRC` + `EGFR/CNA` 双高：高度怀疑 immune-tumor doublet
   * `MBP/PLP1` 极高 + tumor marker 低：可能是 oligodendrocyte 或髓鞘 RNA
   * `SNAP25/SYT1` 少量散在表达：可能是神经元 ambient RNA
   * 用 Scrublet、DoubletFinder、Solo、SoupX、DecontX 辅助处理

4. **空间验证**

   * hypoxia/MES2 应靠近坏死/低氧区域
   * perivascular/ECM 程序应靠近血管或基质区域
   * neuronal/NEU-like tumor 若在浸润边缘且 CNA+，可信度提高
   * 纯正常 neuron/oligodendrocyte 通常位于邻近正常脑或浸润边缘

---

## 6. 最终推荐的 GBM 肿瘤细胞注释字典

| ClassicState | UpdatedState 可选                 | 功能模块可选                                                             | 代表 marker                                                                                    |
| ------------ | ------------------------------- | ------------------------------------------------------------------ | -------------------------------------------------------------------------------------------- |
| AC-like      | AC-like, AC-gliosis, Cilia-like | homeostatic, connectivity, gliosis, cilia                          | `GFAP`, `AQP4`, `ALDOC`, `APOE`, `SLC1A3`, `CHI3L1`, `GAP43`, `DNAH/CFAP`                    |
| OPC-like     | OPC-like, GPC-like              | progenitor, cycling, invasion, oligo-differentiating               | `PDGFRA`, `OLIG1/2`, `SOX10`, `CSPG4`, `PTPRZ1`, `MKI67`, `TOP2A`                            |
| NPC-like     | NPC1-like, NPC2-like, NEU-like  | stem/developmental, neuronal-lineage, cycling, motility            | `SOX2`, `SOX4`, `ASCL1`, `DLL3`, `DCX`, `STMN2`, `NRG1/3`, `NRXN3`, `SYN3`                   |
| MES-like     | MES1-like, MES2/Hypoxia, Stress | inflammatory/gliosis, hypoxia, ECM, invasion, connectivity, stress | `CD44`, `CHI3L1`, `VIM`, `LGALS3`, `ANXA1`, `VEGFA`, `CA9`, `HILPDA`, `TNC`, `FN1`, `HSPA1A` |

---
