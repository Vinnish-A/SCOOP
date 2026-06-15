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

默认工具：FastCore / OmicVerse-backed core engine。当前稳定路径通过 pre-run capability planner 选择后端；未通过验证或环境缺依赖时，唯一 fallback 是 `scanpy_legacy`，即原始 Scanpy normalization/HVG/PCA/kNN/UMAP/Leiden 加 Torch Harmony 路径。

输出：

- H5AD `layers['log1p_norm']`；
- H5AD `obsm`：biology PCA、identity PCA、Harmony PCA、UMAP；
- H5AD `obs`：最终 cluster；
- 外部表：HVG rank、PC covariate association、Leiden sweep、cluster stability。

不做：把所有 resolution/seed 的 cluster 都留在 `obs`；把 UMAP 距离用于注释；对 condition 做 batch correction。

FastCore 额外输出：

- 外部 JSON：`02_core/fastcore/fastcore_manifest.json`、`02_core/fastcore/core_quality.json`；
- 运行日志：记录 selected backend、fallback backend、fallback 是否使用和质量文件路径；
- registry：把 FastCore manifest 和 quality report 注册到 `adata.uns['file_registry']['artifacts']`。

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

默认工具：Scanpy marker；FastDE sparse COSG 可作为快速 marker 辅助；OmicVerse `single.find_markers(method='cosg')` 仅作为可选兼容验证；reference mapping 仅作证据。

输出：

- H5AD `obs`：`cell_type_lvl1/2/3`、`cell_state`、confidence；
- 外部表：cluster markers、annotation evidence、merge/split log。

不做：单独依赖一个 marker database、reference model、GPT 或 CellVote 决定标签。

### 04b_tumor_fastcnvpy

触发条件：样本为肿瘤组织，并且已经完成第一层 `Major Lineage` 注释。

默认逻辑：

- 先用 `cell_type_lvl1` 区分实质性细胞和非实质性细胞；
- 非实质性正常细胞进入 `fastcnv_reference_pool = normal_nonparenchymal`；
- 如果实质性正常细胞能与肿瘤/候选细胞分开，则进入 `normal_parenchymal`；
- 如果正常实质细胞和肿瘤细胞混杂，仍然运行 FastCNVpy，用 pooled nonparenchymal reference 给候选实质细胞提供 CNV 证据；
- FastCNVpy 在合并 H5AD 上构建 pooled reference，但按 `sample_id` 拆分后逐样本计算 CNV，避免单样本 reference 偏移导致误判。

输出：

- H5AD `obs`：`fastcnv_reference_pool`、`fastcnv_cnv_fraction`、`fastcnv_normal_threshold`、`fastcnv_tumor_evidence`；
- 外部表：pooled manifest、pooled cell metadata、pooled chromosome-arm CNV、genomic windows、每样本 FastCNVpy 结果；
- 日志：reference 构成、是否加入正常实质细胞、candidate cell 数、并行参数和 dense/sparse H5AD 模式。

并行策略：默认 `n_jobs=1`，保证峰值内存最低；如果样本数较多且每样本细胞数适中，可以把 `n_jobs` 提到 2-4。每个 worker 只持有一个样本的 selected-gene dense 矩阵，避免把整个合并对象复制到多个进程。大型 H5AD 默认 `--h5ad-mode dense` 追求速度；内存紧张时用 `--h5ad-mode sparse`，但当前基准更慢。

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

默认工具：FastDE pseudobulk DESeq2-like negative-binomial Wald。R DESeq2 和 edgeR quasi-likelihood 保留为 reference validation。

输出：

- 外部表：pseudobulk counts、metadata、FastDE design matrix、FastDE DE table、size factors、dispersions、benchmark/R reference validation table、GSEA；
- H5AD：只存结果索引和必要 summary。

不做：用 Harmony/scVI/integrated expression 做 DE；用 per-cell p value 证明 condition-level 差异；把 marker gene 检验解释成 condition DE；把当前 FastDE 解释成 `lfcShrink` 的完全复刻。

FastDE condition DE 输出列保持 DESeq2 风格：`gene`、`baseMean`、`log2FoldChange`、`lfcSE`、`stat`、`pvalue`、`padj`、`dispersion`，并额外写出 `dispGeneEst`、`dispFit`、`dispMAP`、`dispOutlier`。当前完整 dispersion 基准中，30k genes / 12 pseudobulk samples 下 FastDE 耗时约 `18.55s`，R DESeq2 reference 约 `19.68s`；两者 log2FC Spearman 约 `0.99999`，top100 overlap 为 `1.00`。
