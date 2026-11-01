# 单细胞注释通用流程

`identity_scope` 未声明时仅按 `broad` 记录接受粒度。细层级审定必须明确声明
`identity_scope: fine`；`accepted` 或“已审RNA表型”本身不证明细亚型已解决。
旧表迁移应核对实际审阅证据，不从名称、置信度或单一状态字段推导细层级。

版本：v1.0-refined  
目标：快速、稳健、简洁地完成 scRNA-seq 细胞注释，并把每一步的证据保留到可复核表格中。

---

## 1. 最小输入

```yaml
inputs:
  adata: normalized scRNA-seq object
  clusters: reviewed clustering outputs and within-parent neighbourhoods
  cluster_deg: differential genes for each cluster
  cnmf_programs: batch-handled parent cNMF / FastCNMF programs and usage, required for formal fine review, not broad review
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

若目标是整合图谱，注释主轴必须是当前整合范围的无监督分群与真实基因级
DEG；旧研究和逐来源身份只作参考。逐来源命名后拼接不等于整合后重注释，
也不得用旧标签白名单排除原研究漏识别的竞争谱系。先逐项核验可复用的
QC/分群/CNV/programme checkpoint，再补缺口；参见[多模态准备与复用](../stages/multimodal.md)。
分群及其输入和测量范围核验可靠后，直接在实际表达层计算常规 cluster、
subcluster 和合并后 final-label 基因级 DEG；来源/批次分层是有具体问题时
的可选诊断，不是注释前置条件。核对来源构成、分母及批次混杂的解释限制；
细胞级差异不代替患者/条件层独立重复。
HVG/共同列中缺失的关键marker应回到来源实测基因轴查证，缺列不作为阴性。
marker均值、检测率和特异性排名是描述性证据，不能冒称有真实P/q值的DEG。
原生 `annotate evidence` 导出完整输入基因轴的统计检验；COSG 选项另存辅助
特异性表，不替代 DEG。审阅模板继承登记的 cluster_markers 引用，未引用
分群差异证据的身份决定保留 review_required；提交接口核已登记引用文件
与 hash，仍由分析者审定其生物意义、分母和比较范围。

保留原始带符号检验统计量；例如Wilcoxon/MWU的Z应由U和实际并列秩方差
计算，不能在P下溢为0后经截断P反算Z。后者会使强差异基因大量并列，令top
列表受基因顺序影响。P/q、U、效应量、检测率与排序统计量分开验收；导出
问题可基于已验证U/表达轴补派生统计，不必重拟合有效的整合/图/分群。

一级大类是经过审阅的、有生物意义的细胞类，例如 T/NK，而非数字 cluster
或未经审阅的数字群合并。主进程可依据可靠聚类、常规 DEG 和组织已知
正负 marker 审定大类；大类不以 NMF 为前置。可依据证据合并或拆分数字群，
重审边界并保留竞争及待补证细胞，不要求证明每个细胞纯度。

正式启动细亚群分析前，须有已审定生物父大类和可用的元信息、批次、QC、
测量契约，并完成该父大类的 programme NMF 及适当的真实表达级或模型内部
批次处理。父类审定后，可并行执行 NMF 与准备性类内降维、聚类和常规 DEG；
缺口期间这些工作只能作为补证或探索候选，不得称正式细亚群分析完成。
正式裁定细亚群前还必须集齐第 5 节所有证据；缺一项就记录
`evidence_pending`、主动补证，候选不提交为最终身份。保留完整声明父范围的 UID、来源和竞争上下文；programme
最大 usage 不是身份分类器。探索性数值范围 NMF 不自动满足生物父类要求。

细分稳定性不足时保留每seed与分辨率结果，区分稳定大谱系、边界不确定
和连续状态，不为交付而降低已声明阈值或强命名。先核查实际测量和表示
特征：关键已测基因未进入全图HVG时，可对新整合whole-parent范围重新拟合
可比特征空间；它与原生轴未测基因是两类缺口，不能靠补零或旧标签筛选解决。

继承全图表示后做细分时，先量化实际细型基因在该范围中的信号保留，并结合
真实DEG、原生缺列补证和样本结构判断。保留良好的表示可重建whole新父群的
局部邻域；有具体信息损失才重拟合对应范围，不因一次低ARI默认推翻全部
表示。固定迭代预算的优化不足、连续状态和不可靠身份边界分别记录。

跨范围相同UID的标签须综合实际表达/DEG与竞争证据仲裁，不能按范围名、
旧分流或最大confidence直接覆盖。经审阅的递归子群仅替换它实际细化的父
范围，保留其祖先模型/行号与全部UID；其它独立范围的矛盾仍需审阅。最终
细胞表与合并标签DEG是实际交付，只有命名说明或mapping提案不算提交完成。

```text
Step 0. QC / ambient RNA / doublet review
Step 1. Select final clustering resolution
Step 2. Compute DEG for every cluster
Step 3. Main agent reviews actual clustering, conventional DEG, tissue markers and neighbourhoods to accept meaningful biological broad classes; NMF is not required
Step 4. Freeze the accepted parent UID scope, technical contracts and mixed, conflicting or evidence-pending boundaries
Step 5. Prepare parent batch-handled NMF and within-parent embedding/clustering/DEG, optionally in parallel
Step 6. Start formal fine-subtype analysis only after the parent, technical prerequisites and parent NMF are complete; adjudicate only after all section 5 evidence is complete
Step 7. Review tumor-specific CNV/SNV and context separately from lineage/state when relevant
Step 8. Attach programmes as functional overlays; programme argmax is not an identity classifier
Step 9. Merge redundant reviewed subclusters into canonical final labels when appropriate
Step 10. Recompute DEG per final annotated subtype
Step 11. Export marker tables, annotation evidence table, and quick enrichment report
```

上文cNMF/FastCNMF表示programme任务；实际运行引擎以[programmes](../stages/programs.md)
与receipt为准。原生`sklearn_multistart_nmf`不冒称独立`fastcnmf`。

---

## 3. 决策分支

### 3.1 非肿瘤样本

```text
credible clustering + conventional DEG + tissue markers -> reviewed biological broad class
  -> batch-handled parent NMF and preparatory within-parent clustering/DEG (may be parallel)
  -> completed parent/technical/NMF prerequisites -> formal subtype analysis
  -> all section 5 evidence + main-agent adjudication -> reviewed fine label
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
  -> CNV neutral cells: inspect normal atlas, tumor programmes, patient structure and tissue context; missing evidence is pending, while reviewed conflicting evidence may be unresolved
  -> CNV weak or absent with concordant non-CNV tumor evidence: retain malignant-probable as an option
  -> unresolved multi-source conflict: ambiguous; record alternatives and confidence
```

---

## 4. Broad annotation 规则

大群应是有生物意义且证据支持的一级细胞类（例如 T/NK）；可按实际证据
合并或拆分数值群。下表为谱系线索示例，单个 marker 或数字群不自动成为
最终类；实际大类粒度须由组织背景与竞争证据决定：

```yaml
broad_type:
  immune: PTPRC-positive compartments
  epithelial: EPCAM/KRT/CDH1-positive compartments
  stromal: COL1A1/DCN/LUM/PDGFRA-positive compartments
  endothelial: PECAM1/VWF/KDR/CLDN5-positive compartments
  mural: RGS5/CSPG4/PDGFRB/MYH11-positive compartments
  erythroid_or_platelet: HBB/HBA/PPBP/PF4-positive compartments
# Missing evidence is tracked in workflow progress, not as a broad_type.
```

大群阶段不应该命名到过细亚型，例如：

```text
不要在全体细胞初始分群时直接叫 T_CD8_TEXterm。
不要在未抽取髓系重聚类前直接叫 Mac_TAM_SPP1_TREM2。
不要在未做 CNV 的肿瘤上皮群里直接叫 MaligEpi_LuminalB。
```

---

## 5. Fine annotation 规则

细亚群是有生物意义的群体，而非数字编号。正式分析前先具备已审生物父类、
metadata/batch/QC/measurement 契约及完成批次处理的父类 NMF。候选划分的
类内降维、聚类、UMAP 邻域和常规 DEG 可以作为准备性证据与 NMF 并行生成。
正式裁定前主进程必须亲自审读实际结果，并集齐：

1. 已审父大类及其 UID 范围，真实元信息、批次、QC 和测量可比性及未知项；
2. 父类 NMF 谱、usage、稳定性及真实批次处理和输入尺度；
3. 类内降维、聚类结果与 UMAP 距离/邻域，结合实际表达判断边界；二维 UMAP 距离不能单独决定谱系；
4. 常规 subcluster DEG、正负及竞争谱系 marker、组织知识和可获得的相关背景；
5. 跨样本支持、分母、批次/条件混杂与冲突或缺测的检查；不要求每个亚群存在于所有实验条件，条件特异性可以是真实生物现象；
6. 肿瘤相关时另审 CNV/SNV 或 tumor reference；有合适参考图谱时用作辅助证据。

任一必需项缺失时标记流程状态 `evidence_pending`，主动检索原始来源、组织
marker 和可核验文献，或补跑必要分析；记录哪些证据未获得、不可测或暂不可
执行及原因，继续能做的工作。候选及缺口须显式记录，不能正式命名或提交为
已完成。只有已有足够证据且主进程审阅后仍有竞争、混杂或冲突无法裁定，才
作为身份 `unresolved`；不得用 unresolved 掩盖未检索或未完成的证据。
`evidence_pending` 是工作进度或证据缺口表状态，不强行新增现有 annotation
schema 枚举、不冒充细胞身份。细胞级 marker DEG 用于描述身份，不代替患者/条件效应的独立重复和
相应统计推断；无需为常规聚类 DEG 强制逐来源或批次分层。

cNMF / FastCNMF 的用途：

```text
identity-disambiguation: 去除 cycling/ribosomal/stress 等无意义主轴后，看剩余核心 identity 程序
functional-overlay: 把 IFN、hypoxia、cytotoxic、APC、remodeling 等作为状态后缀
cluster-merge: 若两个小群只有 cNMF 功能强弱不同但历史表型一致，可合并同名
rare-state-detection: 对少量特殊状态给出低/中置信度状态标签，而不是新建“细胞类型”
```

---

## 6. subagent 注释规则

**主进程调度并审核；子进程可执行分析、检索和整合证据并提出候选，最终注释必须由主进程亲自读取实际 DEG、聚类、UMAP 距离/邻域和组织 marker 裁量；细亚群还须亲自读实际父类 NMF 与第 5 节全部证据。不能只接受子代理标签。禁止依赖 SingleR、CellTypist、scANVI 等机械分类或自动标签传播决定最终名称；图谱匹配只作支持材料。逐群记录来源、备选身份、理由和置信度，经结构化验证后提交。治疗响应/生存结局不参与身份命名。CNV 中性本身也不能证明正常，仍需正常身份与组织背景支持。**

子进程可按以下顺序整理候选证据，主进程独立读取实际证据后作最终判断：

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


### 导出审阅表的非空校验

将 DE 统计转入逐群审阅表时，直接对照实际 tested gene axis 和原始
DEG 行。明确区分未测基因、已测但未检出、分母不足未检验，以及未进入
该次请求 panel；不能互相替代。索引为 gene 的表应检查 gene index，
不能误查 DataFrame 列名而把真实已测 marker 记成缺失。

验证须报告实际核对的 group×gene 行数和非空字段数，并独立抽查
正向、负向、极显著和未检验记录的 U、signed Z、P、BH q、均值 FC
及检测率。空集合上的 all/pass 不能作为审阅表已验证的证据。
如果导出逻辑错误使原 DE 证据在审阅中不可见，应归档失效审阅版本并
重审所有受影响身份和状态；原表达与原检验未变时，不为修报告重跑模型。

### 审阅覆盖与身份粒度

结构化 decision 可用 `identity_scope: broad | fine | unresolved`（旧接口默认
fine）。证据缺失先留在工作进度/缺口表并主动补证，不因缺失提交 unresolved；
该状态仅用于已有充分证据但审阅后仍存在混杂或冲突。提交会分别写 `annotation_reviewed`、`broad_annotation_status` 和
`fine_annotation_status`；只审定宽谱系不能自动审定细亚群，也不赋予正常性。
被拒绝的新 decision 必须撤销对应细胞旧的 accepted 资格，保留旧文字供追溯。
明确区分所有记录已审、宽谱系已定、细亚群已定、QC 已测与参考可用。

分群抽样的证据可能掩盖小的异质残群。对最终未定 UID 子集单独核查竞争谱系
和抽样覆盖；有系统性竞争信号时回到真实原生特征、无监督分群和全基因 DEG。
不能把漏查、panel 未请求、HVG 未含某谱系基因写成数据天然无法辨认。
