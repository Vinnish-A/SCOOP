# scSP Agent SOP Project

这是一个面向 Agent 执行的单细胞 / 空间转录组 SOP 工程模板。它不是“工具大全”，而是一套可执行、可审计、可复用的项目结构。

本项目的核心标准是：**快速、稳健、简洁**。这里的“简洁”不是把说明文字省掉，而是把概念和交付物压缩到必要范围：每个模块只承担一类责任，每类任务只有一个默认方法，大型结果以表格或 sidecar 文件保存，H5AD 只保存核心状态和结果索引。

## 概念层面的简洁

项目只保留七个分析模块：

1. `01_qc`：逐样本 QC、Scrublet 双细胞、ambient 标记或修正。
2. `02_core`：normalization、HVG、PCA、Torch Harmony、kNN、UMAP、Leiden sweep。
3. `03_programs`：快速 consensus NMF，必要时用 OmicVerse cNMF 验证。
4. `04_annotation`：marker、NMF、kNN、reference、ontology 融合注释。
5. `05_spatial`：多分辨率空间结构、RCTD-py 反卷积、空间 graph。
6. `06_ccc`：FastCCC 初筛、CellPhoneDB/LIANA 多聚体验证。
7. `07_de`：pseudobulk DE 与下游统计。

这七个模块覆盖常规项目的大部分需求，但避免引入 trajectory、velocity、CNV、foundation model zoo、drug response、full OmicVerse lazy pipeline 等非默认分析。

## H5AD 克制原则

H5AD 是主状态对象，不是数据库。它保存：

- 每个 cell / spot / bin 的核心状态字段；
- 当前使用的少数 embedding 和 graph；
- gene flags；
- 外部结果文件的 registry。

FastCCC 结果、DE 长表、marker 表、NMF gene weights、Leiden sweep、空间几何、多边形、transcript coordinates、模型 checkpoint 和完整审计日志都写在 `runs/<run_id>/` 中，不塞进 `adata.uns`。

## OmicVerse 的定位

OmicVerse 是复用层，而不是 SOP 主控层。项目只复用它提供的基础设施和成熟 wrapper：

- I/O：`omicverse.io.read_h5ad`、`read_10x_h5`、`read_10x_mtx`、`read_visium_hd`、`read_xenium`、`save`。
- GPU/CPU 转换：`omicverse.pp.anndata_to_GPU`、`anndata_to_CPU`。
- marker wrapper：`omicverse.single.find_markers`，尤其是 COSG 用作快速 marker 辅助。
- cNMF wrapper：`omicverse.single.cNMF`，作为 NMF 稳健性验证层。
- CCC 验证：`omicverse.single.run_cellphonedb_v5`、`run_liana`。
- 轻量报告：`omicverse.single.generate_scRNA_report`，作为报告草稿。

不把 `omicverse.single.lazy`、CellVote、GPTCelltype、trajectory zoo、velocity zoo、drug response、CNV zoo 等纳入默认 SOP。

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

# 7. 运行空间 RCTD-py 反卷积（空间项目）
python scripts/04_spatial_rctd.py \
  --config runs/2026-06-12_001/config/run.yaml

# 8. 运行 CCC 分析
python scripts/05_ccc_fastccc.py \
  --config runs/2026-06-12_001/config/run.yaml

# 9. 瘦身 H5AD 并生成 manifest
python scripts/06_prune_and_manifest.py \
  --config runs/2026-06-12_001/config/run.yaml
```

NMF 程序发现默认使用 `programs.method=fastcnmf`，`max_iter=50`；OmicVerse cNMF 仅作为显式验证 fallback。

## 关键配置文件

- `project.yaml`：项目级路径、样本字段、存储策略。
- `configs/default_run.yaml`：默认 run 配置。
- `configs/h5ad_schema.yaml`：H5AD 中允许保留的字段。
- `configs/omicverse_reuse_policy.yaml`：哪些 OmicVerse 设施可用，哪些禁用。
- `docs/`：SOP、接口说明、参数调优规则。
- `src/scsp_agent_sop/`：可复用 Python 模块。
- `scripts/`：Agent 可直接调用的命令行脚本。
