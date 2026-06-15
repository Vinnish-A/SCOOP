
## 2. 决策分级：先身份，再生态位，再功能

| 层级 | 决策问题                               | 强标志物/排除项                                                                                                                             | 输出                                                                    |
| -- | ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------- |
| L0 | 是否为成纤维细胞谱系？                        | 阳性：`COL1A1/COL1A2/COL3A1/DCN/LUM/COL6A1/DPT/PDGFRA`；排除：`PTPRC` 免疫、`EPCAM/KRT` 上皮、`PECAM1/VWF` 内皮、`RGS5/CSPG4/MYH11` 高表达 pericyte/SMC | `FB` 或排除                                                              |
| L1 | 正常/驻留 fibroblast 还是 activated CAF？ | normal-like：`PI16/DPP4/CD34/MFAP5/COL15A1/DPT`；activated CAF：`FAP/PDPN/POSTN/MMP11/COL11A1/CTHRC1/LRRC15`                            | `resident-FB` 或 `CAF`                                                 |
| L2 | 属于哪个生态位？                           | 结合空间位置、邻近细胞、配体-受体、组织结构标志                                                                                                             | `N:vascular / tumor-border / matrix / TLS / crypt / adipose / immune` |
| L3 | 主功能模块是什么？                          | ECM、contractile、inflammatory、apCAF、IFN、reticular、vascular、hypoxic/glycolytic、cycling 等模块分数                                           | `F:...`                                                               |
| L4 | 癌种特异细化                             | CRC 保留 crypt/WNT 轴；乳腺癌保留 perivascular/adipose/CAF-S1 等；PDAC 保留 myCAF/iCAF/apCAF/LRRC15 轴                                             | `CRC-...` / `BRCA-...`                                                |
| L5 | 证据等级                               | A：转录组 + 蛋白/空间；B：多 marker + pathway；C：单 marker 或弱证据                                                                                   | `E:A/B/C`                                                             |

尤其要注意：**FAP、ACTA2、PDPN 不能单独定义 CAF subtype**。乳腺癌单细胞研究指出缺乏一个共同 CAF marker；Krishnamurty 等也显示 `Fap` 和 `Acta2` 可见于 fibroblasts 和 pericytes，而 `LRRC15` 对 TGFβ-driven CAF 更具限制性。([Nature][3])

## 3. CAF 功能—生态位词典 v1.0

| 模块                       | 建议名称                              | 强标志物                                                                      | 主要功能                               | 生态位判断                                |
| ------------------------ | --------------------------------- | ------------------------------------------------------------------------- | ---------------------------------- | ------------------------------------ |
| Matrix CAF               | `F:ECM-remodeling`                | `MMP11/CDH11/POSTN/COL1A1/2/3/COL11A1/FN1/TNC`                            | ECM 沉积、胶原重塑、侵袭前沿                   | stroma/invasive front                |
| Contractile myCAF        | `F:contractile`                   | `ACTA2/TAGLN/MYL9/TPM1/2/CNN1/CTGF`                                       | 收缩、基质硬化、TGFβ 反应                    | tumor-proximal 或 desmoplastic stroma |
| LRRC15+ myCAF            | `F:TGFβ-immunosuppressive matrix` | `LRRC15/POSTN/THBS2/COL11A1/CTHRC1`                                       | CD8 T 细胞抑制、ICB 抵抗、病理性基质 setpoint   | tumor-proximal matrix                |
| Inflammatory iCAF        | `F:cytokine/chemokine`            | `IL6/LIF/CXCL1/2/8/CXCL12/CCL2/PLA2G2A/CFD/C3/CD34`                       | 炎症、趋化、JAK/STAT、免疫调节                | 免疫浸润区、血管附近                           |
| Complement/secretory CAF | `F:complement-inflammatory`       | `C3/C7/CFB/CFD/CFH/CFI`                                                   | 补体调节、炎症放大                          | immune-rich 或 tumor-adjacent         |
| Antigen-presenting CAF   | `F:apCAF`                         | `CD74/HLA-DRA/HLA-DRB1/HLA-DPA1/HLA-DPB1/CIITA`；必须 `PTPRC−/LST1−/CD68−`   | MHC-II 抗原呈递/免疫调节                   | immune niche；需排除 myeloid doublet     |
| IFN-response CAF         | `F:IFN-response`                  | `IDO1/ISG15/MX1/OAS1/IFIT1/2/3/CXCL9/10/11`                               | IFN 反应、免疫激活/抑制混合状态                 | tumor-stroma border 或 inflamed tumor |
| Vascular CAF             | `F:vascular support`              | `MCAM/CD146/NOTCH3/COL18A1/ACTA2low-int`，pericyte 应 `RGS5/CSPG4/PDGFRB` 高 | 血管生成、血管稳定                          | PECAM1/VWF+ 内皮邻近                     |
| Reticular/TLS CAF        | `F:reticular/TLS organization`    | `CCL19/CCL21/CXCL13/IL7/PDPN/LTBR`                                        | TLS 形成、T/B 细胞组织化、T cell attraction | TLS/B cell follicle 周围               |
| Tumor-border/hypoxic CAF | `F:tumor-like/hypoxic-glycolytic` | `MME/CD10/NT5E/CD73/ENO1/NDRG1/CA9/HSPH1/PDPN`                            | 肿瘤边界适应、缺氧/糖酵解、肿瘤接触                 | tumor-stroma border                  |
| Cycling/stress CAF       | `S:cycling` 或 `S:stress`          | `MKI67/TOP2A/UBE2C`；`HSPA1A/HSPH1/DNAJB1`                                 | 增殖或应激状态                            | 作为状态后缀，不单独当 lineage                  |

Cords 等的 Table 1 已经把 matrix、inflammatory、vascular、tumour-like、IFN-response、apCAF、reticular-like、dividing CAF 的 marker、预期功能和空间分布并列列出，例如 matrix CAF 为 `MMP11/CDH11/POSTN/collagens`，inflammatory CAF 为 `PLA2G2A/CFD/C3/CD34/CXCL12/CD248`，reticular-like CAF 为 `CCL21/CCL19` 且位于 TLS 周围。([Nature][4]) PDAC 研究进一步支持 iCAF 的 `IL6/IL8/CXCL1/CXCL2/CXCL12/CCL2` 模块、myCAF/iCAF 的可塑性，以及 apCAF 的 MHC-II/CD74 模块。([PMC][1])

## 4. 癌种场景中的具体落地

### 结直肠癌 CRC

CRC 不应只注释为 “myCAF/iCAF”。至少保留四类层级：

1. **`CRC-CAF | N:matrix/invasive | F:ECM-remodeling`**
   对应早期 CRC 单细胞研究中的 CAF-A，强 marker：`MMP2/DCN/COL1A2`，功能为 ECM 重塑。([PMC][1])

2. **`CRC-CAF | N:vascular-contractile | F:contractile/myCAF`**
   对应 CAF-B，强 marker：`ACTA2/TAGLN/PDGFA`；但若 `RGS5/CSPG4/PDGFRB/MYH11` 高，应优先判为 pericyte/SMC，而不是 CAF。([PMC][1])

3. **`CRC-FB/CAF | N:crypt-bottom or crypt-top | F:epithelial niche`**
   结肠本身有强生态位轴：crypt-bottom fibroblasts 为 `PDGFRAlow`，分泌 `WNT2/WNT2B/RSPO3/GREM1` 支持 stem cell；crypt-top fibroblasts 为 `PDGFRAhigh`，分泌 `WNT5A` 和 BMP ligands 促进上皮分化。这个轴在 CRC 中不能被粗暴并入普通 iCAF 或 myCAF。([PLOS][5])

4. **`CRC-CAF | N:tumor-border/hypoxic | F:hypoxic-glycolytic`**
   若有 `CA9/NDRG1/ENO1/MME/NT5E`，且空间上靠近肿瘤边界，应作为 tumor-border/hypoxic CAF，而不是仅称 “activated CAF”。

### 乳腺癌 BRCA

乳腺癌应保留 **perivascular、matrix、developmental/adipose、immune-suppressive** 等生态位和功能。Bartoschek 等在乳腺癌模型中定义三类空间和功能上不同的 CAF，并将其与 perivascular niche、mammary fat pad 和 transformed epithelium 等来源/位置联系起来。([Nature][3])

推荐最小标签集：

* **`BRCA-CAF | N:perivascular | F:vascular support | markers:NID2/MCAM/NOTCH3/COL18A1`**
* **`BRCA-CAF | N:collagen-rich stroma | F:ECM-remodeling | markers:PDGFRA/MMP11/POSTN/COL1A1`**
* **`BRCA-CAF | N:immune | F:CXCL12-mediated immunosuppression | markers:FAP/PDPN/CXCL12/IL6/PDGFRB`**
* **`BRCA-CAF | N:TLS-adjacent | F:reticular/T-cell attraction | markers:CCL19/CCL21/PDPN`**
* **`BRCA-CAF | N:tumor-border | F:tumor-like/hypoxic | markers:MME/NT5E/NDRG1/ENO1/CA9`**

乳腺癌中的 CAF-S1 是一个很好的例子：它不是简单的 “iCAF”，而是具备 CXCL12 介导 CD4+CD25+ T 细胞募集、OX40L/PD-L2/JAM2 介导保留、B7H3/CD73/DPP4 促进 Treg 生存和分化的免疫抑制功能模块。([PubMed][6]) 后续 CAF-S1 单细胞研究又拆出 8 个 cluster，其中 ECM-myCAF 和 TGFβ-myCAF 与免疫治疗原发耐受相关，提示 “CAF-S1” 也应再分解为功能子模块，而不是作为终点标签。([PubMed][7])

### PDAC 及其他实体瘤

PDAC 可以作为功能验证最充分的模板：`myCAF` 通常为 ACTA2/TAGLN/CTGF 高、靠近肿瘤细胞；`iCAF` 为 IL6/LIF/CXCLs 高、相对远离肿瘤细胞；`apCAF` 为 `CD74/MHC-II` 高；`LRRC15+ myCAF` 是 TGFβ-driven、免疫抑制、限制 anti-PD-L1 反应的强功能 subtype。([PMC][1]) 对 NSCLC、HNSCC、卵巢癌、胃癌、胰腺癌等，建议使用同一套功能轴，但必须保留癌种/组织后缀与生态位证据。最近 pan-cancer 空间多组学研究在 10 个癌种、超过 1400 万细胞中发现 4 类保守空间 CAF subtype，并显示其邻域、免疫表型和临床结局相关，说明“空间生态位”应成为 CAF 注释的一级属性，而非后期补充。([PubMed][8])

## 5. 证据等级与冲突规则

**E:A，高可信。** 至少 2–3 个强 marker 同向；module score 高；有空间/IHC/IMC/RNAscope 或蛋白验证；能说明邻近细胞和功能。例如 `CCL19/CCL21/PDPN` 且位于 TLS 周围，可标为 `reticular/TLS CAF`。

**E:B，中可信。** 多 marker + pathway 支持，但无空间验证。例如 `MMP11/POSTN/COL1A1/COL11A1` 高，可标为 `matrix CAF`，但生态位写成 `N:matrix-inferred`。

**E:C，低可信。** 单 marker、dropout 明显、疑似 doublet 或不同模块冲突。输出为 `CAF-like | F:unresolved`，不要强行命名。

冲突处理：
`MKI67/TOP2A` 高只加 `S:cycling` 后缀；`HSPH1/HSPA1A` 高只加 `S:stress/heat-shock` 后缀；`CD74/HLA-DRA` 高但同时 `PTPRC/LST1/CD68` 高，应先判定免疫 doublet 或 myeloid contamination；`ACTA2/TAGLN` 高但 `RGS5/CSPG4/MYH11` 高，应优先考虑 pericyte/SMC。
