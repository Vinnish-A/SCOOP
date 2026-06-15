# 03. OmicVerse 复用策略

OmicVerse 当前是一个统一的 Python 多组学平台，覆盖 bulk、single-cell、spatial、visualization、model-based analysis 和 AI-assisted workflows。它的价值在于提供很多 AnnData-native 的设施和 wrapper，而不是替代本 SOP 的决策逻辑。

## 推荐复用

### I/O 和数据接入

可复用：

- `omicverse.io.read_h5ad`
- `omicverse.io.read_10x_h5`
- `omicverse.io.read_10x_mtx`
- `omicverse.io.read_visium_hd`
- `omicverse.io.read_visium_hd_bin`
- `omicverse.io.read_visium_hd_seg`
- `omicverse.io.read_xenium`
- `omicverse.io.save`

这些函数可以减少不同平台输入格式的胶水代码，但读入后仍要进入本 SOP 的 schema validation。

### GPU/CPU 转换

可复用：

- `omicverse.pp.anndata_to_GPU`
- `omicverse.pp.anndata_to_CPU`

用途：大规模对象在 GPU/CPU 间切换时减少重复代码。但这不是算法选择，不改变 QC、HVG、batch correction、clustering 的规则。

### Marker 辅助

可复用：

- `omicverse.single.find_markers`
- `omicverse.single.get_markers`

`find_markers` 支持 COSG、wilcoxon、logreg 等方法。COSG 可作为快速 marker 发现，适合大对象快速浏览；正式 evidence table 仍需要写出为外部表并进入 annotation fusion。

### NMF 验证

可复用：

- `omicverse.single.cNMF`

本 SOP 默认用 FastCNMF profile。只有 programme 稳定性不足、关键结论依赖 programme、或需要 publication-grade 验证时才调用 OmicVerse cNMF wrapper。

### CCC 验证

可复用：

- `omicverse.single.run_cellphonedb_v5`
- `omicverse.single.run_liana`
- `omicverse.single.format_liana_results`
- OmicVerse CCC 可视化函数

FastCCC 是默认初筛。CellPhoneDB v5 / LIANA 只用于多聚体、机制关键 LRI 或排序稳健性验证。

### 报告草稿

可复用：

- `omicverse.single.generate_scRNA_report`

它可以生成 MultiQC-style 报告草稿，但不能替代本 SOP 的审计日志和模块报告。

## 不建议作为默认

不纳入默认：

- `omicverse.single.lazy`
- CellVote / GPTCelltype 作为最终注释器
- trajectory zoo
- velocity zoo
- drug response
- CNV zoo
- perturbation zoo
- foundation model zoo
- OmicVerse full preprocessing pipeline

这些功能不是坏工具，但它们会扩大 SOP 概念范围，使常规项目交付复杂化。只有当具体生物问题要求时才作为扩展分析运行。
