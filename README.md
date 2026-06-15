# SCOOP: Single Cell Omics Operating Protocol

SCOOP（Single Cell Omics Operating Protocol）是一个面向单细胞 / 空间转录组分析的可执行、可审计、可复用操作协议工程。它不是“工具大全”，而是把常见分析任务整理成一套 Agent 和分析人员都能稳定执行的项目结构。

本项目的核心标准是：**快速、稳健、简洁**。这里的“简洁”不是把说明文字省掉，而是把概念和交付物压缩到必要范围：每个模块只承担一类责任，每类任务只有一个默认方法，大型结果以表格或 sidecar 文件保存，H5AD 只保存核心状态和结果索引。

## 架构层

SCOOP 当前按四层组织：

1. **SOP Workflow Layer**：run 目录、模块顺序、H5AD 状态、file registry、decision log。
2. **Fast Compute Layer**：确定性计算引擎，包括 FastDE、FastCNMF、FastCNVpy、FastCopyKAT；SCOOP 通过 `src/scoop_fast/` 提供统一 contract 和 registry。
3. **Evidence & Skill Layer**：`markers/` 和未来 `skills/` 中的 marker、state、tumor、naming 规则，只作为证据来源。
4. **Annotation Decision Layer**：`src/scsp_agent_sop/annotation_decision/` 负责 evidence bundle、结构化 decision schema、validator 和 committer。

Agent/AI 层不能直接改 H5AD 或自由发明标签。它只能产生结构化 annotation decision；只有通过确定性 validator 的 decision 才会由 committer 写入 H5AD。

## 概念层面的简洁

项目主流程保留七个分析模块，其中 `04_annotation` 内部拆成 evidence、tumor CNV gate、decision draft、validated commit 四个步骤：

1. `01_qc`：逐样本 QC、Scrublet 双细胞、ambient 标记或修正。
2. `02_core`：normalization、HVG、PCA、Torch Harmony、kNN、UMAP、Leiden sweep。
3. `03_programs`：FastCNMF gene programme discovery，必要时用 OmicVerse cNMF 显式验证。
4. `04_annotation`：marker/NMF/reference/CNV evidence、annotation decision validation、commit。
5. `05_spatial`：多分辨率空间结构、RCTD-py 反卷积、空间 graph。
6. `06_ccc`：FastCCC 初筛、CellPhoneDB/LIANA 多聚体验证。
7. `07_de`：FastDE pseudobulk DESeq2-like condition DE 与下游统计。

这七个模块覆盖常规项目的大部分需求，但避免引入 trajectory、velocity、CNV、foundation model zoo、drug response、full OmicVerse lazy pipeline 等非默认分析。

## H5AD 克制原则

H5AD 是主状态对象，不是数据库。它保存：

- 每个 cell / spot / bin 的核心状态字段；
- 当前使用的少数 embedding 和 graph；
- gene flags；
- 外部结果文件的 registry。

FastCCC 结果、DE 长表、marker 表、NMF gene weights、Leiden sweep、空间几何、多边形、transcript coordinates、模型 checkpoint、annotation decision audit 和完整审计日志都写在 `runs/<run_id>/` 中，不塞进 `adata.uns`。

## Fast 引擎

SCOOP 暴露并复用这些确定性 Fast engine：

- `fastde markers`：sparse COSG-like marker scoring，不 densify 大矩阵。
- `fastde deseq2`：pseudobulk DESeq2-like negative-binomial Wald test。
- `fastcnmf`：FastCNMF programme discovery，默认 `n_iter=20`、`max_iter=50`。
- `fastcnvpy`：pooled-reference tumor CNV evidence。
- `fastcopykat`：CNV prediction compatibility path。

`src/scoop_fast/registry.py` 只记录这些 engine 的 contract、输入输出 schema、CLI 和质量门；它不重写算法。

## OmicVerse 的定位

OmicVerse 是可选复用层，不是 SOP 主控层。当前默认 marker 快速路径是 FastDE sparse COSG，OmicVerse `single.find_markers(method='cosg')` 仅作为兼容验证；`single.cNMF` 只作为 NMF 稳健性验证 fallback。SCOOP 不把 `omicverse.single.lazy`、CellVote、GPTCelltype、trajectory zoo、velocity zoo、drug response、CNV zoo 等纳入默认 SOP。

## 快速开始

```bash
# 1. 安装项目本身
pip install -e .

# 2. 复制 run 模板
cp -r runs/template runs/2026-06-12_001

# 3. 编辑配置
vim runs/2026-06-12_001/config/run.yaml

# 4. 运行 QC
python scripts/01_qc_scrublet.py \
  --config runs/2026-06-12_001/config/run.yaml

# 5. 运行核心分析
python scripts/02_core_analysis.py \
  --config runs/2026-06-12_001/config/run.yaml

# 6. 运行 NMF 程序发现
python scripts/03_fast_consensus_nmf.py \
  --config runs/2026-06-12_001/config/run.yaml

# 7. 导出 annotation evidence，不直接赋最终标签
python scripts/04_annotation_markers.py \
  --config runs/2026-06-12_001/config/run.yaml

# 8. 肿瘤项目可运行 pooled-reference FastCNVpy evidence gate
python scripts/04b_tumor_fastcnvpy.py \
  --config runs/2026-06-12_001/config/run.yaml \
  --gene-metadata data/external/references/gene_metadata.tsv

# 9. 生成 evidence bundle 和 annotation decision template
python scripts/04c_annotation_decide.py \
  --config runs/2026-06-12_001/config/run.yaml

# 10. 人工或 Agent 编辑 structured decision JSON 后，验证并提交标签
python scripts/04d_annotation_commit.py \
  --config runs/2026-06-12_001/config/run.yaml \
  --decisions runs/2026-06-12_001/04_annotation/decisions/annotation_decision_template.json

# 11. 运行空间 RCTD-py 反卷积（空间项目）
python scripts/05_spatial_rctd.py \
  --config runs/2026-06-12_001/config/run.yaml

# 12. 运行 CCC 分析
python scripts/06_ccc_fastccc.py \
  --config runs/2026-06-12_001/config/run.yaml

# 13. 准备 pseudobulk 并运行 condition DE
python scripts/08_prepare_pseudobulk.py \
  --config runs/2026-06-12_001/config/run.yaml
fastde deseq2 runs/2026-06-12_001/07_de/pseudobulk/<cell_type> condition ctrl test

# 14. 瘦身 H5AD 并生成 manifest
python scripts/07_prune_and_manifest.py \
  --config runs/2026-06-12_001/config/run.yaml
```

NMF 程序发现默认使用 `programs.method=fastcnmf`，`max_iter=50`；OmicVerse cNMF 仅作为显式验证 fallback。Annotation commit 只提交 validator 接受的 decision；需要 review 的 cluster 保持 `annotation_status=review_required`。

## 关键配置文件

- `project.yaml`：项目级路径、样本字段、存储策略。
- `configs/default_run.yaml`：默认 run 配置。
- `configs/annotation_decision_schema.yaml`：annotation decision 允许字段、confidence、tumor CNV gate 和 program sanitizer 规则。
- `configs/h5ad_schema.yaml`：H5AD 中允许保留的字段。
- `configs/omicverse_reuse_policy.yaml`：哪些 OmicVerse 设施可用，哪些禁用。
- `docs/`：SOP、接口说明、参数调优规则和架构说明。
- `markers/`：人工整理的 marker/state/tumor/naming 参考资料。
- `src/scoop_fast/`：Fast engine contract layer。
- `src/scsp_agent_sop/`：SCOOP 的可复用 Python 模块。
- `scripts/`：Agent 可直接调用的命令行脚本。
