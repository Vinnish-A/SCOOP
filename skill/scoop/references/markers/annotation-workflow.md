# 单细胞注释通用流程

版本：v1.0-refined  
目标：快速、稳健、简洁地完成 scRNA-seq 细胞注释，并把每一步的证据保留到可复核表格中。

---

## 1. 最小输入

```yaml
inputs:
  adata: normalized scRNA-seq object
  clusters: final cluster labels after best clustering resolution is selected
  cluster_deg: differential genes for each cluster
  cnmf_programs: cNMF / FastCNMF programs and per-cell usage matrix
  cnv_result: FastCNVpy result for tumor-relevant compartments, optional for non-tumor
  project_context:
    species: human | mouse | other
    tissue: <tissue>
    disease: tumor | normal | inflamed | other
    cancer_type: BRCA | CRC | GBM | other
    sample_info: patient, region, treatment, pathology if available
```

---

## 2. 总体流程

```text
Step 0. QC / ambient RNA / doublet review
Step 1. Select final clustering resolution
Step 2. Compute DEG for every cluster
Step 3. Broad cell type annotation by marker + exclusion
Step 4. Tumor-specific FastCNVpy branch if tumor sample
Step 5. Recluster each broad compartment when fine annotation is needed
Step 6. Compute DEG for every subcluster
Step 7. Run cNMF / FastCNMF within each broad compartment
Step 8. Annotation subagent maps marker evidence to historical phenotypes
Step 9. cNMF programs are attached as functional overlays, not primary identities
Step 10. Merge redundant subclusters into canonical final labels when appropriate
Step 11. Recompute DEG per final annotated subtype
Step 12. Export marker tables, annotation evidence table, and quick enrichment report
```

---

## 3. 决策分支

### 3.1 非肿瘤样本

```text
cluster DEG + marker evidence
  -> broad cell type
  -> compartment-level reclustering
  -> subcluster DEG
  -> cNMF/FastCNMF functional programs
  -> skill historical phenotype mapping
  -> final label
```

### 3.2 肿瘤样本，独立审阅恶性证据

```text
broad annotation
  -> identify epithelial / parenchymal / glial / other candidate malignant compartments
  -> CNV/SNV, normal reference and patient/tissue-context review
  -> malignant_status assignment separate from differentiation/state
  -> cancer-type references support candidate review; malignant subtype names require supported malignancy
  -> normal tissue skill for normal-reference/context-supported parenchymal cells; CNV-neutral alone does not route cells to normal
```

#### 情况 A：正常实质与肿瘤细胞已经分开

FastCNVpy 提供辅助证据：

```text
tumor-like programme + concordant CNV/context -> reviewed malignant support
coherent normal lineage + matched normal-reference + tissue-context support -> normal-high; CNV neutrality alone is insufficient
```

#### 情况 B：正常实质与肿瘤细胞没有分开

联合证据区分来源：

```text
same marker-defined parenchymal/epithelial cluster
  -> CNV abnormal cells: supporting evidence; review reference, programme and context
  -> CNV neutral cells: inspect normal atlas, tumor programmes, patient structure and tissue context; may be normal, malignant-probable or unresolved
  -> CNV weak or absent with concordant non-CNV tumor evidence: retain malignant-probable as an option
  -> unresolved multi-source conflict: ambiguous; record alternatives and confidence
```

---

## 4. Broad annotation 规则

大群注释只需要做到：

```yaml
broad_type:
  immune: PTPRC-positive compartments
  epithelial: EPCAM/KRT/CDH1-positive compartments
  stromal: COL1A1/DCN/LUM/PDGFRA-positive compartments
  endothelial: PECAM1/VWF/KDR/CLDN5-positive compartments
  mural: RGS5/CSPG4/PDGFRB/MYH11-positive compartments
  erythroid_or_platelet: HBB/HBA/PPBP/PF4-positive compartments
  unknown: evidence insufficient
```

大群阶段不应该命名到过细亚型，例如：

```text
不要在全体细胞初始分群时直接叫 T_CD8_TEXterm。
不要在未抽取髓系重聚类前直接叫 Mac_TAM_SPP1_TREM2。
不要在未做 CNV 的肿瘤上皮群里直接叫 MaligEpi_LuminalB。
```

---

## 5. Fine annotation 规则

精细亚群注释需要同时查看：

1. subcluster DEG；
2. 历史表型 marker；
3. 与相邻谱系的排除 marker；
4. cNMF / FastCNMF 程序；
5. 样本、组织、病理、治疗、空间位点；
6. 参考图谱映射结果；
7. 在肿瘤细胞中必须查看 CNV/SNV 或 tumor reference。

cNMF / FastCNMF 的用途：

```text
identity-disambiguation: 去除 cycling/ribosomal/stress 等无意义主轴后，看剩余核心 identity 程序
functional-overlay: 把 IFN、hypoxia、cytotoxic、APC、remodeling 等作为状态后缀
cluster-merge: 若两个小群只有 cNMF 功能强弱不同但历史表型一致，可合并同名
rare-state-detection: 对少量特殊状态给出低/中置信度状态标签，而不是新建“细胞类型”
```

---

## 6. subagent 注释规则

**最终注释由分析者或推理模型结合组织匹配图谱、可核验原始文献、实际marker/排除marker、相邻群和跨样本证据综合判断。禁止依赖SingleR、CellTypist、scANVI等机械分类或自动标签传播来决定最终细胞名称；图谱匹配只作支持材料。逐群记录来源、备选身份、理由和置信度，经结构化验证后提交。治疗响应/生存结局不参与身份命名。CNV中性本身也不能证明正常，仍需正常身份与组织背景支持。**

subagent 读取每个 cluster/subcluster 的证据后，按如下顺序执行：

```text
1. Check broad lineage gate.
2. Reject incompatible lineages using exclusion markers.
3. If tumor-relevant, check malignancy_status from FastCNVpy/SNV/reference.
4. Search matching historical phenotypes in the relevant marker skill.
5. If a canonical historical phenotype is supported, use it.
6. If phenotype is uncertain but real cell identity is clear, use <real phenotype + function>.
7. If two small clusters map to the same real phenotype and only differ by generic cNMF programs, merge naming.
8. Record confidence and do_not_call labels.
```

---

## 7. 置信度

```yaml
high:
  required:
    - broad lineage marker set consistent
    - at least two to three core phenotype markers consistent
    - key exclusion markers absent or low
    - tissue / cancer / CNV / reference context consistent
medium:
  required:
    - broad lineage reliable
    - phenotype marker set partially consistent
    - some uncertainty in exclusion, CNV, or reference mapping
low:
  triggered_by:
    - single-marker call
    - ambient RNA or doublet risk
    - CNV weak or absent where required
    - conflicting historical phenotypes
    - very small cluster or low UMI
```

---

## 8. 最终导出

每轮注释结束后必须导出：

```text
cluster_markers.tsv/parquet
subcluster_markers.tsv/parquet
final_label_markers.tsv/parquet
annotation_evidence.tsv
cnmf_program_summary.tsv
fastcnvpy_summary.tsv, if tumor
quick_enrichment_by_final_label.tsv
```

`final_label_markers` 应在最终命名合并之后重新计算，不能只沿用原始 cluster DEG。
