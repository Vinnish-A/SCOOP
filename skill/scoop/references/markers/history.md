# markers skill 重构记录

版本：v1.0-refined

---

## 1. 结构变化

- 新增 `annotation_workflow.md`：把项目级注释流程、FastCNVpy 分支、cNMF/FastCNMF 的位置写成 SOP。
- 新增 `naming_strategy.md`：明确“属 + 种差 + 状态后缀 + 证据等级”的通用命名策略。
- 新增 `major_cell_types.md`：提供大群注释 marker 词典和排除项。
- 新增 `tumor-general.md`：把肿瘤样本中 malignant / normal 判定抽象成通用规则。
- 新增 `FINAL_NAMING_REFERENCE.md`：给出最终可查阅命名表。

---

## 2. 内容变化

- T 细胞改为 2025 modular nomenclature 风格：历史表型 + marker + modular 属性，不再把单一状态过度本体化。
- Myeloid 删除 M1/M2 作为最终标签，改为 monocyte / macrophage / DC / TAM state 的分层注释。
- Fibroblast / CAF 保留真实 CAF 表型和癌种生态位，但把 ECM、contractile、IFN、cycling 等作为状态或功能 overlay。
- BRCA / CRC / GBM 肿瘤细胞按癌种单独处理；所有 malignant subtype 都要求先有 malignant_status。
- CNV 用于 malignant/normal 与克隆证据，不直接解释功能状态。
- cNMF / FastCNMF 用于功能连续谱和 identity disambiguation，不作为 skill 内部硬阈值打分器。

---

## 3. 删除或降级的内容

- 删除基于单基因打分直接命名的倾向。
- 删除或降级“功能就是细胞类型”的表述。
- `Cycling`、`RibosomalHigh`、`Stress`、`Hypoxiahi`、`IFNhi` 等统一降为 overlay。
- `CMS1-4` 不再作为 CRC 单细胞 epithelial label；只作为 sample/region pseudobulk 解释层。
- `Normal-like` 不再作为 BRCA 肿瘤主标签；默认进入 normal/ambient/doublet/CNV 复核。
- GBM 的 `Cycling` 不作为第五经典状态；作为 AC/OPC/NPC/MES 上的 overlay。

---

## 4. subagent 预期行为

```text
读入 marker evidence + DEG + cNMF programs + FastCNVpy + project context
  -> 查找对应 markers skill
  -> 先判 identity 和 malignant_status
  -> 再映射历史 phenotype
  -> 最后挂载 functional overlay
  -> 输出 final_label、confidence、uncertainty、do_not_call
```
