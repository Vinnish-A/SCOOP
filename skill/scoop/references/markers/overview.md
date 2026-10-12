# 单细胞注释 markers skill 总说明

版本：v1.0-refined  
适用范围：当前项目的单细胞转录组学注释，不覆盖空间、蛋白、多组学验证的完整实验规范。  
核心定位：`markers/` 下的文件不是自动分类器，也不是基因打分脚本，而是给注释 subagent 使用的“证据词典 + 决策流程”。

---

## 1. 本次重构目标

本目录只保留三类信息：

1. **历史上真实可发表的细胞表型或命名体系**  
   例如 `T_CD8_TEMRA`、`T_CD8_TEXp`、`Mac_TAM_SPP1_TREM2`、`CRC_MalignantEpi_iCMS2`、`GBM_Malignant_AC-like`、`BRCA_MalignantEpi_LuminalA`。

2. **这些表型在单细胞转录组中可见的 marker / signature 证据**  
   marker 应作为“共同证据集合”使用，不应单基因一票定名。

3. **从粗注释到精细注释的决策顺序**  
   每个 skill 文件都包含 `decision_process`、`required_evidence`、`do_not_call` 和最终命名建议。

功能解释、细胞周期、核糖体、缺氧、IFN、应激、代谢等连续谱，原则上由 cNMF / FastCNMF 提供。markers skill 只把 cNMF 程序映射为**可叠加状态后缀**，不再在 skill 内部设计打分算法。

---

## 2. 文件结构

| 文件 | 作用 |
|---|---|
| `annotation_workflow.md` | 项目级注释流程：DEG、cNMF、FastCNVpy、subagent 决策、最终导出 |
| `naming_strategy.md` | 通用命名策略：属 + 种差 + 状态后缀 + 证据等级 |
| `major_cell_types.md` | 大群注释 marker 与排除项 |
| `T.md` | T 细胞历史表型、marker 与 2025 modular nomenclature 决策规则 |
| `myeloid.md` | 单核/巨噬/DC/粒细胞/肥大细胞注释规则 |
| `fibroblast.md` | fibroblast / CAF 注释规则 |
| `tumor-general.md` | 肿瘤样本通用 CNV 与 malignant/normal 判定规则 |
| `tumor-BRCA.md` | 乳腺癌恶性上皮与正常乳腺实质上皮注释 |
| `tumor-CRC.md` | 结直肠癌恶性上皮与正常肠上皮注释 |
| `tumor-GBM.md` | GBM 恶性胶质瘤细胞与正常脑实质细胞注释 |
| `FINAL_NAMING_REFERENCE.md` | 可直接查阅的最终命名参考表 |

---

## 3. 全局原则

### 3.1 先身份，后功能

- `cell identity`：由 marker、排除项、参考映射、CNV/SNV、组织背景共同判断。
- `functional state`：由 cNMF / FastCNMF、DEG、功能富集和上下文解释；不应覆盖身份。
- `technical state`：cycling、stress、ribosomal、mitochondrial、ambient RNA、doublet 风险只作为附加说明。

### 3.2 大群注释与亚群注释不同

- 大群注释可主要依赖 marker 与排除项。
- 亚群精细注释必须联合：
  - 小群 DEG；
  - cNMF / FastCNMF 程序；
  - 参考文献中的历史表型；
  - 样本、组织、病理、CNV/SNV 证据；
  - 必要时把相近小群合并为同一规范大群或同一历史表型。

### 3.3 肿瘤组织必须先解决 malignant / normal

在肿瘤样本中，实质细胞与肿瘤细胞可能分不开。处理顺序为：

1. 先完成大群注释。
2. 对实质/上皮/胶质等可疑肿瘤来源细胞运行 FastCNVpy。
3. 若正常实质与肿瘤细胞已经分开：FastCNVpy 用作辅助确认。
4. 若没有分开：联合 FastCNVpy 连续证据、正常图谱、肿瘤程序和患者/组织背景判断恶性；CNV低/中性不能自动划入正常。
5. 只有 malignant 细胞进入癌种专属 tumor skill。
6. 正常实质细胞按对应组织的正常细胞谱系命名。

### 3.4 不确定就是输出不确定

证据不足时输出：

```text
Unknown_<broad_type>
Ambiguous_<broad_type>_<reason>
<cell_type>_like_lowConfidence
<cell_type>_<historical_phenotype>_mixed
```

禁止把“只有一个 marker 高”的 cluster 强行命名为精细细胞亚型。

---

## 4. 推荐输出 schema

```yaml
annotation_record:
  cluster_id: <cluster>
  broad_type: <major lineage>
  parent_type: <canonical parent cell type>
  historical_phenotype: <published phenotype if supported>
  functional_overlay: <cNMF-derived state, optional>
  malignancy_status: malignant | normal | probable_malignant | ambiguous | not_applicable
  final_label: <canonical final name>
  display_label: <short label for figures>
  positive_evidence: [marker genes, protein markers, CNV/SNV, reference labels]
  exclusion_evidence: [markers excluding adjacent lineages]
  cnmf_evidence: [program ids and top genes, if used]
  confidence: high | medium | low
  uncertainty: <remaining alternatives or reason>
  do_not_call: [labels explicitly rejected]
```
