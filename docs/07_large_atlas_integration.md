# 07. 大型单细胞图谱整合：scVI、可信注释与 scPoli

状态：SOP 设计文档；质控与预处理复用 SCOOP 现有模块，本文只规定大型图谱的整合、注释和验证门槛。

## 1. 目标与边界

本文用于将多个 donor、sample、study 或测序批次整合成可查询的大型单细胞图谱。目标是：

- 用 scVI 学习无监督、批次感知的潜在表示；
- 在 `X_scVI` 中重建细胞结构并重新聚类；
- 在 scVI 阶段内完成 marker、外部参考、跨样本复现性、cNMF 功能和人工证据整合，产出图谱标签；
- 生成一份可供用户逐群审核的 scVI annotation report；
- 只让高置信度标签参与 scPoli 的监督原型学习；
- 将中低置信度细胞保留为 unlabeled，并用留出的完整 donor 或 study 做验证。

本文不重新定义 QC 阈值，不用 integrated expression 做 DE，也不把 UMAP 几何当作注释证据。QC、doublet、ambient RNA 和 suspicious-cell 处理沿用 [`qc`](../templates/run/qc/README.md) 与 [`单细胞注释通用流程`](../markers/annotation_workflow.md)。condition-level DE 仍使用 raw counts 派生的 pseudobulk，而不是 `X_scVI` 或 scPoli latent space。

## 2. 本流程的硬约束

### 2.1 scVI 模型不直接生成标签，但 scVI 阶段必须产出标签

scVI 输出批次感知的潜在表示，例如：

```python
adata.obsm["X_scVI"] = scvi_model.get_latent_representation()
```

`X_scVI` 中的 Leiden 编号只是 cluster，不是有生物学意义的细胞类型。但在本 SOP 中，“scVI 阶段”不是只运行模型，而是一个完整交付单元：

```text
scVI model
  -> X_scVI
  -> graph / clustering
  -> marker / reference / distribution / neighbour / cNMF review
  -> hierarchical annotation
  -> confidence
  -> scVI annotation report
```

因此 scVI 阶段结束时必须产出 `cell_type_lvl1/2/3`、`cell_state`、置信度和审核报告。准确含义是：**标签由基于 scVI 空间的完整 annotation 流程产出，而不是由 scVI 模型头直接预测。**

### 2.2 弃用 scANVI 和旧低质量注释

本流程不再把 scANVI 作为大型图谱的默认整合或标签传播步骤。旧注释若来源不清、粒度不一致、缺少 marker 支持或已知质量较低，不得作为 scPoli 的监督标签。

旧标签只允许保留在独立字段中用于审计和差异比较，例如：

```text
obs["legacy_cell_type"]
obs["legacy_annotation_source"]
```

它们不得覆盖新标签，不得参与 `labeled_indices`，也不得作为新标签正确性的唯一证据。若某一旧标签经过当前 marker、参考和跨 donor 复核，可把这些新证据登记为一次新的 annotation decision；这不等于直接继承旧标签。

### 2.3 scPoli 前必须人工确认

scPoli 的 prototype loss 会使用输入标签约束潜在空间。错误标签可能形成错误 prototype，并将错误传播到更多细胞。因此，scPoli 不是自动修复低质量标签的工具。

**在人工或项目负责人明确接受 scVI 阶段标签质量前，流程必须停在 annotation gate，不得训练 scPoli。**

## 3. 总体流程

```text
各数据集原始 counts
    ↓
复用 SCOOP qc：逐样本 QC / doublet / ambient RNA / suspicious flags
    ↓
统一 gene、cell 和 metadata schema，保留 counts
    ↓
按 donor / sample / study 预先划分 train、validation、held-out query
    ↓
scVI 无监督整合
    ↓
在 X_scVI 上构图、Leiden、UMAP 和稳定性分析
    ↓
marker + 负 marker + 外部参考 + 跨 donor/study 一致性
    + 图邻近群比较 + cNMF identity/function programmes
    ↓
分层注释，并赋予 high / medium / low confidence
    ↓
生成 scVI annotation report
    ↓
强制 annotation gate：用户根据报告确认标签质量
    ├── reject / revise → 返回聚类或注释
    └── approve
          ↓
仅 high-confidence cells 作为 scPoli labeled cells
medium / low / ambiguous cells 保持 unlabeled
          ↓
scPoli atlas training
          ↓
完整 held-out donor / study query mapping 与开放集检查
```

## 4. 输入与数据冻结

### 4.1 最小输入

每个细胞至少需要：

```yaml
required:
  counts: adata.layers["counts"]
  sample_id: technical sample identifier
  donor_id: biological replicate identifier
  study_id: source study or cohort
  condition: disease / treatment / control, not a default correction key
  species: human | mouse | other
recommended:
  library_id: library or lane
  chemistry: assay chemistry
  tissue: tissue or anatomical site
  dataset_id: immutable source dataset identifier
  legacy_cell_type: audit only
```

`cell_id` 必须在全图谱内唯一，推荐写成 `<dataset_id>:<original_barcode>`。基因需要统一到同一物种和同一标识体系；重复 symbol、跨物种同源映射、缺失 gene ID 和不同 genome build 必须在训练前记录。

### 4.2 先冻结验证单元

在查看 scPoli 结果前按完整 biological unit 划分：

- training：用于 scVI 结构发现、人工注释和 scPoli 训练；
- validation：用于调 prototype 距离阈值和检查已知类型映射；
- held-out query：完整 donor 或 study，仅在最终验证时作为 unlabeled query 使用。

禁止随机按 cell 切分作为主要验证，因为同一 donor 的近邻细胞会同时进入训练和验证，产生过度乐观的结果。稀有类型若无法完整留出，应在报告中明确标记，而不是悄悄改成 cell-level random split。

## 5. 复用 QC 与预处理

每个来源数据集先执行 SCOOP `qc`，并复用已有输出：

- `obs` 中的 QC metrics、`doublet_score`、`doublet_call` 和 `final_use`；
- `qc_thresholds_by_sample.tsv`；
- `scrublet_summary_by_sample.tsv`；
- `scrublet_scores.parquet`；
- ambient RNA、low-quality 和 suspicious-cell flags。

合并大型图谱时只做以下额外检查：

1. 所有数据集的 `final_use` 语义一致；
2. raw counts 未被 log-normalized matrix 替代；
3. QC 删除比例按 study、donor、condition 汇总，避免某一 biological group 被系统性移除；
4. doublet 和 ambient RNA flag 在合并后仍可追溯到原 sample；
5. 高 mt、低 genes、stress、ribosomal 和 cell-cycle 是否在 `X_scVI` 中形成技术主轴。

大型整合不应为了“批次混合更好”而放回已确认的 doublet 或 QC failure。对于身份清楚但带有 stress、IFN、cycling 等状态的细胞，沿用 SCOOP 原规则：优先保留细胞并把程序记录为 `cell_state`，不要把它们命名为新 cell type。

## 6. scVI 无监督整合

### 6.1 模型职责

scVI 负责学习表达分布和潜在表示；邻居图、Leiden、UMAP 和注释都在模型训练之后单独完成。最小示意代码如下：

```python
import scanpy as sc
import scvi

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="sample_id",
)

scvi_model = scvi.model.SCVI(adata)
scvi_model.train()

adata.obsm["X_scVI"] = scvi_model.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI", key_added="scvi")
sc.tl.leiden(adata, neighbors_key="scvi", key_added="leiden_scvi")
sc.tl.umap(adata, neighbors_key="scvi")
```

代码只表达流程，不固定版本参数。运行时必须记录 `scvi-tools`、PyTorch、CUDA、数据 schema、随机种子、训练参数和模型 checkpoint。

### 6.2 batch key 规则

默认只把技术来源作为 batch/covariate 候选，例如 `sample_id`、`library_id`、`chemistry` 或明确的 processing batch。不得默认校正 `condition`、`disease`、`treatment`、`response` 或目标 lineage。

当 study 与 disease 完全混淆时，模型无法区分 study effect 和真实 disease effect。此时不得用“UMAP 混合良好”宣称整合成功；应保留混淆说明，并通过同 condition 的跨 study 数据、held-out study 和 marker 保真度做敏感性验证。

### 6.3 大型图谱训练检查

至少检查：

- training/validation loss 是否收敛且无数值异常；
- 每个 study、donor 和 major lineage 是否在 latent space 中有合理覆盖；
- kNN 是否被单一 study、library size、mt/ribosomal/stress 程序主导；
- canonical identity marker 是否仍保真；
- 已知 biology 是否被过度混合，技术重复是否仍完全分离；
- 不同随机种子或合理超参数下的邻域和 cluster 是否稳定。

UMAP 只用于展示。整合质量判断必须结合 latent kNN、marker、sample composition、cluster stability 和 held-out mapping。

## 7. 在 `X_scVI` 上建立可信标签

### 7.1 分层注释

第一轮只标 major lineage，例如 T、B、myeloid、epithelial、stromal 和 endothelial。第二轮在每个 major lineage 内重新构图或聚类，再判断细亚型。肿瘤候选实质细胞继续复用 SCOOP FastCNVpy 分支，不能仅凭 scVI/scPoli latent space 判断 malignant status。

每个 cluster 的命名至少综合：

- canonical marker 和 negative marker；
- cluster differential expression；
- 至少一个适配组织与物种的外部参考；关键类型建议使用多个独立参考；
- cluster 是否跨 donor、sample 和 study 重复出现；
- 相邻 cluster 的可区分性及不同 seed/resolution 下的稳定性；
- cNMF/FastCNMF 的 identity programme 与功能 programme；
- doublet、low quality、ambient RNA、stress、cell cycle 和 disease-state 混杂；
- 稀有类型是否有足够细胞与独立 donor 支持。

参考映射是证据之一，不是最终裁决。参考之间冲突时必须降低 confidence 并记录 conflict。

### 7.2 标签与置信度字段

推荐写入：

```text
obs["cell_type_lvl1"]
obs["cell_type_lvl2"]
obs["cell_type_lvl3"]
obs["cell_state"]
obs["scvi_curated_cell_type"]  # 当前细胞可接受的最细粒度标签
obs["annotation_confidence"]   # high | medium | low
obs["annotation_status"]       # labeled | unlabeled | ambiguous | excluded
obs["annotation_source"]       # scvi_cluster_curated
obs["annotation_version"]
```

置信度和 scPoli 输入规则：

| 证据状态 | annotation | scPoli 用法 |
| --- | --- | --- |
| marker 明确、负 marker 无冲突、跨 donor/study 稳定 | `high` | 作为 labeled cell |
| 大类明确、亚型证据不足 | `medium` | 在细亚型模型中保持 unlabeled；粗粒度需另作 high-confidence decision |
| marker 冲突、过渡态、单 donor 稀有群或技术风险 | `low` | 保持 unlabeled |
| doublet、确认的 QC failure | `excluded` | 不进入训练 |

不要把所有 cluster 强制命名。不要把 `unknown`、`ambiguous`、`low_quality` 当成普通 cell type 送入监督学习，否则模型会为异质混合群学习一个 prototype。

### 7.3 scVI annotation report

scVI 阶段必须同时输出机器可读表格和人类可读报告。推荐主文件为：

```text
reports/scvi_annotation_report.html
```

如需要归档，可同时导出 PDF；PDF 不替代可筛选、可跳转的 HTML 和 sidecar tables。报告是用户批准 scPoli annotation gate 的主要依据，不是可选展示图。

#### A. 报告首页

首页给出：

- 输入细胞、QC 后细胞、study、donor、sample 和 condition 数量；
- scVI 的 batch/covariate 设置、模型版本、checkpoint 和训练状态；
- major lineage 与亚群数量；
- high/medium/low/unlabeled/excluded 的细胞数和比例；
- 拟进入 scPoli `labeled_indices` 的细胞数、标签数和 donor/study 覆盖；
- 未解决的 annotation conflicts、单 donor 亚群和 stop reasons；
- 全图谱 UMAP，分别按 cluster、最终标签、study、donor、condition、confidence 和主要 QC flag 着色。

#### B. 亚群分布

每个亚群都要同时报告 cell count 和组成比例，不能只给 UMAP：

- 在 study、donor、sample、condition、tissue 和 chemistry 中的细胞数；
- `P(metadata_group | subtype)`：该亚群由哪些来源构成；
- `P(subtype | metadata_group)`：该亚群在每个来源中的相对丰度；
- 出现该亚群的 donor/study 数量及最小支持阈值；
- 是否由单一 donor、study、condition 或 batch 主导；
- 低细胞数来源需标记，不能把抽样噪声解释为生物学缺失。

推荐同时提供 stacked bar、dot plot 和长表。condition 差异这里只作描述性展示；正式丰度推断仍应使用 donor/sample-level abundance analysis。

#### C. 邻近群

每个亚群列出最邻近的 3-5 个 cluster/label，并解释为什么保留拆分或建议合并。邻近关系必须来自 `X_scVI` kNN graph，例如 cluster 间归一化边权、shared-neighbour strength 或 latent centroid distance；不得根据二维 UMAP 上看起来接近来定义。

每一对邻近群至少展示：

- source、target、双向 kNN edge weight 和 rank；
- latent centroid distance，作为辅助而非唯一指标；
- 共同 marker、区分 marker 和相互 negative marker；
- donor/study composition 差异；
- cNMF programme usage 差异；
- `keep_split | merge | uncertain` 结论及理由。

若两个亚群只有 stress、IFN、cycling、ribosomal 或 hypoxia 等连续功能 programme 强弱不同，而 identity marker 和分布一致，优先使用同一 cell type 加不同 `cell_state`，不要自动创造两个 cell types。

#### D. cNMF/FastCNMF 功能

复用 SCOOP `programs` 的 FastCNMF 结果；必要时在 major lineage 内重跑，而不是直接用全图谱中由 lineage 构成决定的 programme。每个 programme 先分类为：

```text
identity
functional_state
technical_or_nuisance
mixed_or_unresolved
```

每个亚群报告：

- dominant programmes 及 usage 的 median、IQR 和高 usage 细胞比例；
- 每个 programme 的 top genes 和简短生物学解释；
- programme 在 donor、study、condition 中是否复现；
- 它支持的是 cell identity、cell state，还是技术干扰；
- 相邻亚群间 programme usage 的差异；
- programme 与最终标签、`cell_state` 和 confidence 的关系。

cNMF programme 不能单独决定 cell type。由 cell cycle、stress、ribosomal、mitochondrial、dissociation 或 ambient RNA 主导的 programme 必须标为 state/nuisance，并在 identity 注释时降权。

#### E. 每群 annotation card

报告为每个 cluster/subcluster 生成一张可独立审核的 annotation card，字段至少包括：

```yaml
cluster_id: <id>
proposed_label: <cell_type_lvl1/2/3>
cell_state: <state or null>
confidence: high | medium | low
n_cells: <int>
n_donors: <int>
n_studies: <int>
distribution_summary: <text + linked table>
nearest_clusters: <3-5 graph neighbours>
positive_markers: <genes>
negative_markers: <genes and conflicts>
reference_support: <references and concordance>
cnmf_identity_programmes: <programmes>
cnmf_functional_programmes: <programmes>
qc_and_doublet_risk: <summary>
alternative_labels: <labels or null>
decision: keep | merge | split | unlabeled | exclude
reason: <auditable rationale>
review_required: true | false
```

annotation card 应链接原始 marker、programme、distribution、neighbour 和 reference tables，避免只显示截断后的 top genes。

#### F. 报告 sidecar tables

报告至少由以下表格支撑：

```text
tables/scvi_subgroup_distribution.parquet
tables/scvi_cluster_neighbours.parquet
tables/scvi_cluster_markers.parquet
tables/scvi_annotation_cards.tsv
tables/scvi_annotation_conflicts.tsv
tables/cnmf_program_classification.tsv
tables/cnmf_usage_by_subgroup.parquet
tables/cnmf_usage_by_donor_study.parquet
tables/reference_mapping_concordance.tsv
```

`scvi_cluster_neighbours.parquet` 至少含 `source_cluster`、`target_cluster`、`edge_weight`、`edge_rank`、`centroid_distance` 和 `mutual_neighbour`。`scvi_subgroup_distribution.parquet` 使用长表，至少含 `subgroup`、`metadata_key`、`metadata_value`、`n_cells`、`fraction_within_subgroup`、`fraction_within_metadata_value`、`n_donors` 和 `n_studies`。

## 8. scPoli 前强制 annotation gate

### 8.1 必须提交的 review bundle

运行 scPoli 前应导出：

```text
annotation/atlas/
  scvi_model_manifest.json
  scvi_cluster_summary.tsv
  scvi_cluster_markers.parquet
  scvi_annotation_report.html
  scvi_subgroup_distribution.parquet
  scvi_cluster_neighbours.parquet
  cnmf_program_classification.tsv
  cnmf_usage_by_subgroup.parquet
  annotation_decisions.tsv
  annotation_conflicts.tsv
  annotation_confidence_by_cluster.tsv
  annotation_confidence_by_donor.tsv
  label_by_study_donor.tsv
  reference_mapping_concordance.tsv
  figures/scvi_umap_by_label_study_donor_qc.pdf
  scpoli_annotation_gate.yaml
```

review bundle 必须允许审核者逐群看到亚群分布、graph neighbours、marker/negative marker、cNMF identity/功能、各 donor/study 构成、reference 一致性、QC 风险和拟用于 supervision 的细胞数。HTML 报告及其全部 sidecar tables 共同计算 `review_bundle_hash`。

### 8.2 go/no-go 清单

用户或项目负责人需要逐项确认：

- [ ] scVI 阶段已通过完整 annotation 产出标签，且没有把 Leiden 编号当作 cell type；
- [ ] 已逐群阅读 scVI annotation report 中的分布、邻近群、marker、cNMF 功能和冲突；
- [ ] scANVI 和未经复核的 legacy labels 未参与当前监督标签生成；
- [ ] major lineage 与细亚型粒度一致，命名遵循 SCOOP schema；
- [ ] high-confidence label 有 marker、negative marker 和跨样本证据；
- [ ] 中低置信度、冲突、过渡态和可疑细胞未进入 `labeled_indices`；
- [ ] 每个受监督 label 的 donor/study 支持已检查，单 donor 类型已单独说明；
- [ ] validation 与 held-out query 按完整 donor/study 冻结；
- [ ] 肿瘤标签引用了 CNV/SNV/病理或适用的肿瘤参考证据；
- [ ] 标签没有被 stress、cell cycle、ribosomal 或 ambient RNA 主轴定义；
- [ ] 审核者接受当前标签错误可能被 scPoli prototype loss 放大。

gate 文件建议最少包含：

```yaml
schema_version: "1.0"
run_id: <run_id>
label_version: <annotation_version>
status: pending  # pending | approved | rejected
approved_by: null
approved_at: null
review_bundle_hash: <sha256>
notes: null
```

只有 `status: approved`、审核人和时间完整、bundle hash 与当前结果一致时才能继续。标签、筛选规则或 `labeled_indices` 发生变化后，旧批准自动失效，必须重新审核。

### 8.3 必须回退的情况

出现以下任一情况时不得运行 scPoli：

- high-confidence 标签仍主要继承 legacy annotation；
- 同一 label 的核心 marker 在大部分 donor 中缺失；
- 一个 label 只由 batch、study、condition 或通用状态程序定义；
- reference、marker 和 lineage gate 存在未解决冲突；
- 拟监督的 label 粒度混杂，例如一部分是 major lineage、一部分是极细亚型；
- held-out donor/study 尚未冻结；
- review bundle 不完整或 gate 未明确批准。

## 9. scPoli 半监督训练

### 9.1 只监督 high-confidence cells

概念示意：

```python
import numpy as np

supervised = (
    adata.obs["annotation_confidence"].eq("high")
    & adata.obs["annotation_status"].eq("labeled")
    & adata.obs["final_use"].astype(bool)
)
labeled_indices = np.flatnonzero(supervised.to_numpy())

# 下面的构造参数必须按项目锁定的 scarches 版本核对。
scpoli_model = scPoli(
    adata=adata,
    condition_keys=["sample_id"],
    cell_type_keys=["scvi_curated_cell_type"],
    labeled_indices=labeled_indices,
)
```

`medium` 不参与监督。如果细亚型不确定但 major lineage 可靠，应另行记录该细胞在粗粒度层级的 high-confidence decision，训练独立的 coarse atlas；不能把不确定细亚型伪装成 high confidence。若项目混合了不同注释粒度，优先分别训练 coarse atlas 和 lineage-specific fine atlas。

### 9.2 原型与类别检查

训练前后都应报告：

- 每个 label 的细胞数、donor 数和 study 数；
- 是否存在单 donor prototype 或严重类别不平衡；
- 同一标签是否需要多个 prototype；
- 标签内部 latent 异质性和相邻标签分离度；
- unlabeled cell 到最近 labeled prototype 的距离分布；
- 不同 seed 下 prototype 和映射结果的稳定性。

对过大的常见类型可做 donor-aware 下采样，避免 prototype loss 被大类支配；对过小类型不得仅通过重复采样制造“证据”。采样规则和随机种子必须写入 manifest。

## 10. 验证与避免循环论证

不能用“scPoli 成功预测了训练标签”证明标签正确。标签已参与 prototype learning，这只能说明模型复现了监督信号。

### 10.1 donor/study-level 验证

将完整 validation 或 held-out donor/study 作为 unlabeled query 映射，检查：

- predicted label 与 canonical/negative marker 是否一致；
- 已知类型在 donor 和 study 间的 sensitivity、precision 与 confusion；
- query 是否被迫映射到训练集中不存在的错误类别；
- 疾病特异状态是否被误判为新 identity；
- 最近 prototype 距离、预测不确定性和 donor/study drift；
- 稀有类型在独立 biological replicate 中是否复现。

prototype distance 没有通用 cutoff。阈值应在 validation donors 上按 cell type、类别规模和错误成本校准，并在 held-out query 上冻结使用。阈值之外的细胞应返回 unknown/unassigned，而不是强制最近邻标签。

### 10.2 最低验收条件

一个 atlas release 至少满足：

1. 所有监督标签都有可追溯的新证据和 confidence；
2. scPoli 训练使用的确切 `labeled_indices` 可复现；
3. annotation gate 已批准且 hash 匹配；
4. 至少一个完整 donor 或 study 未参与监督训练；
5. held-out 结果同时报告准确性、unknown rate 和 prototype distance；
6. 主要 label 在多个 donor 中保留 marker coherence；
7. known failure、单 donor 类型、参考冲突和低置信度区域写入 release notes。

## 11. 输出与审计

### 11.1 H5AD 中保留的最小状态

```text
layers["counts"]
obsm["X_scVI"]
obsm["X_scPoli"]                 # 训练后，如版本接口提供
obs["leiden_scvi"]
obs["cell_type_lvl1/2/3"]
obs["cell_state"]
obs["scvi_curated_cell_type"]
obs["annotation_confidence"]
obs["annotation_status"]
obs["atlas_split"]               # train | validation | held_out
obs["scpoli_supervised"]
obs["scpoli_predicted_label"]
obs["scpoli_prototype_distance"]
uns["file_registry"]              # 指向外部模型、表格、图和 manifest
```

字段名需按实际 scarches 版本和 SCOOP schema 对齐。长表、模型 checkpoint、全量距离矩阵和训练历史不写入 H5AD。

### 11.2 外部交付物

```text
models/scvi/<model_version>/
models/scpoli/<model_version>/
tables/scvi_cluster_markers.parquet
tables/annotation_decisions.tsv
tables/scpoli_labeled_indices.parquet
tables/prototype_summary.tsv
tables/heldout_mapping_metrics.tsv
tables/prototype_distance_calibration.tsv
reports/scvi_annotation_report.html
reports/atlas_release_notes.md
logs/decision_log.jsonl
```

所有关键文件登记 path、schema version、software version 和 hash。atlas release 必须同时记录 source dataset 版本、gene universe、QC 版本、annotation version、scVI model version、scPoli model version 和 held-out split。

## 12. 推荐执行顺序

1. 各 study/sample 复用 SCOOP `qc`，统一 counts 与 metadata schema；
2. 冻结 train、validation 和 held-out donor/study；
3. 训练 scVI，导出 `X_scVI`、模型 manifest 和 checkpoint；
4. 在 `X_scVI` 上做稳定性驱动的 Leiden 聚类；
5. 先注释 major lineage，再在各 major lineage 内细化亚型；
6. 整合亚群分布、graph neighbours、marker/reference 和 cNMF identity/功能，产出正式标签与 high/medium/low confidence；
7. 生成 scVI annotation report 和 review bundle，等待用户明确批准 annotation gate；
8. 仅把 high-confidence cells 放进 scPoli `labeled_indices`；
9. 将 medium/low/ambiguous cells 作为 unlabeled，与高置信度细胞共同训练或映射；
10. 用完整 held-out donor/study 验证 marker、标签、unknown detection 和 prototype distance；
11. 验收后冻结模型、标签、阈值、hash 和 release notes。

## 13. 参考资料

- [scVI user guide](https://docs.scvi-tools.org/en/1.3.3/user_guide/models/scvi.html)
- [scPoli model interface](https://github.com/theislab/scarches/blob/master/scarches/models/scpoli/scpoli_model.py)
- [scPoli population-level integration paper](https://www.nature.com/articles/s41592-023-02035-2)
- [scPoli surgery and label-transfer pipeline](https://docs.scarches.org/en/stable/scpoli_surgery_pipeline.html)
