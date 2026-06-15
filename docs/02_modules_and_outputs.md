# 02. 模块与交付物

## 01_qc

输入：`adata.layers['counts']`、`obs['sample_id']`、gene symbols。

默认工具：Scanpy QC metrics，逐样本 robust MAD/quantile threshold，`scanpy.pp.scrublet`。

输出：

- H5AD `obs`：QC 指标、QC flag、doublet score、doublet call；
- 外部表：`qc_thresholds_by_sample.tsv`、`scrublet_scores.parquet`；
- 图：QC violin/scatter、Scrublet histograms；
- 日志：每个 sample 的阈值和删除比例。

不做：全局 `mt > 10%` 删除；默认删除 ribosomal-high 或 proliferation-high；把 Scrublet 分数当成唯一 doublet 证据。

## 02_core

输入：QC 后 H5AD。

默认工具：Scanpy normalization/HVG/PCA/kNN/UMAP/Leiden，Torch Harmony。

输出：

- H5AD `layers['log1p_norm']`；
- H5AD `obsm`：biology PCA、identity PCA、Harmony PCA、UMAP；
- H5AD `obs`：最终 cluster；
- 外部表：HVG rank、PC covariate association、Leiden sweep、cluster stability。

不做：把所有 resolution/seed 的 cluster 都留在 `obs`；把 UMAP 距离用于注释；对 condition 做 batch correction。

## 03_programs

输入：log-normalized expression 或 counts 派生的非负矩阵。

默认工具：FastCNMF profile，使用 exact coordinate-descent NMF、K sweep、20 个 replicate seeds 和 `max_iter=50`。OmicVerse `single.cNMF` 只作为稳健性验证。

输出：

- H5AD `obsm['X_nmf_usage']`，如果程序数和 cell 数可控；
- 外部表：NMF gene weights、K sweep、programme summary；
- H5AD `obs`：dominant programme 和 programme entropy 可选。

不做：默认运行所有 NMF zoo；把单次 NMF 结果当作稳定 gene program。

## 04_annotation

输入：final clusters、marker table、NMF usage、kNN graph、reference mapping。

默认工具：Scanpy marker；OmicVerse `single.find_markers(method='cosg')` 可作为快速 marker 辅助；reference mapping 仅作证据。

输出：

- H5AD `obs`：`cell_type_lvl1/2/3`、`cell_state`、confidence；
- 外部表：cluster markers、annotation evidence、merge/split log。

不做：单独依赖一个 marker database、reference model、GPT 或 CellVote 决定标签。

## 05_spatial

输入：spatial H5AD view、reference H5AD、坐标和空间 metadata。

默认工具：RCTD-py for deconvolution，Squidpy/spatial graph，BANKSY for domain 可选。

输出：

- H5AD `obs`：spatial unit metadata、deconvolution summary、domain label；
- H5AD `obsm['spatial']`；
- 外部表：RCTD weights、mode sensitivity、spatial domain table；
- SpatialData/Zarr 或 GeoParquet：图像、polygon、transcript points。

默认规则：低分辨率 spot/ROI 用 RCTD-py `full`；如果每个单位预期只有少数主导 cell types，再运行 `multi` 作为敏感性分析。

## 06_ccc

输入：annotated H5AD、cell-type group、LR resource。

默认工具：FastCCC。复杂 ligand/receptor 用 CellPhoneDB v5 或 LIANA 验证，可通过 OmicVerse wrapper 调用。

输出：

- 外部表：FastCCC all/significant/differential，complex-sensitive validation；
- H5AD `uns`：方法摘要、路径、hash、显著结果数量；
- 可选 H5AD `obs`：cell-type-level sender/receiver summary。

不做：把 FastCCC 全结果写进 `adata.uns`；把 CCC 当作机制证明；忽略多聚体 ligand/receptor 的 subunit 表达。

## 07_de

输入：raw counts、cell type/state label、donor/sample/condition metadata。

默认工具：pseudobulk edgeR quasi-likelihood。

输出：

- 外部表：pseudobulk counts、metadata、design matrix、DE table、GSEA；
- H5AD：只存结果索引和必要 summary。

不做：用 Harmony/scVI/integrated expression 做 DE；用 per-cell p value 证明 condition-level 差异。
