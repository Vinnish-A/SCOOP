# 00. 设计原则

## 目标

SCOOP（Single Cell Omics Operating Protocol）为 Agent 和分析人员执行 single cell omics 项目提供一套稳定路线，覆盖单细胞和空间组学中的常规交付任务。它不是完整的生信百科，而是把常规项目压缩成少数确定模块。每个模块要清楚说明：输入是什么、默认工具是什么、输出写到哪里、什么情况下调参、什么情况下需要人工 review。

## 真正的简洁

文本简短不等于项目简洁。真正的简洁来自概念边界清楚：

- QC 只决定哪些观测是可靠、可疑或失败；
- core analysis 只产生稳定 embedding 和 cluster；
- programs 只解释连续状态和干扰程序；
- annotation 只融合证据，不创造未经证实的细胞类型；
- spatial 只处理空间单位、图和反卷积；
- CCC 只产生候选通信，不证明机制；
- DE 只做统计检验，不用 integrated expression。

## 默认方法唯一

每类任务只设一个默认方法：

- doublet：`scanpy.pp.scrublet`，逐样本运行；
- batch correction：Harmony 2.0；
- NMF：fast consensus NMF；
- spatial deconvolution：RCTD-py；
- low-resolution RCTD-py mode：`full`；
- CCC：FastCCC；
- complex-sensitive CCC validation：CellPhoneDB v5 或 LIANA；
- DE：FastDE pseudobulk DESeq2-like NB Wald；R DESeq2/edgeR only as reference validation.

Fallback 不是并列默认。只有默认失败或触发具体风险时才使用 fallback，并写入 decision log。

## H5AD 不是数据库

H5AD 只保存当前状态。长表和大对象写出到 run 目录。Agent 每次写结果时必须判断：

1. 这个结果是否是每 cell / spot / bin 的核心状态？
2. 它是否下游经常用于过滤、分组或展示？
3. 它是否足够小且稳定？
4. 它是否必须和 AnnData 一起传输？

如果答案不是明确的“是”，就写外部表格，并把路径、hash 和 schema 登记到 `adata.uns['file_registry']`。
