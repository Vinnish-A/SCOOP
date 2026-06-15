# 05. 参数调整规则

参数调整的原则是：先判断问题来自数据、模型还是生物学，再修改最少的参数。

## QC

不要用全局固定阈值。每个 sample 用自己的分布。只有以下情况触发人工 review：

- 一个 sample 被 QC fail 超过 30%；
- ribosomal-high 超过 20% 且影响 PCA；
- high mt 和 low genes 不共现，但被单指标判为 fail；
- 某 condition 因 QC 几乎消失；
- 空间 FOV/tile 与 QC failure 完全重合。

## Batch correction

默认 Harmony 2.0。只校正技术 key，例如 `sample_id`、`library_id`、`batch_id`、`chemistry`。不得校正 `condition`、`disease`、`treatment` 或 lineage。

如果 batch 与 condition 高度混淆，不做 correction，而在下游 DE 或 stratified analysis 中处理。

## Leiden resolution

不要只跑一个 resolution。默认 sweep `0.2–2.0`，5 个 seed。选最低足够 resolution：

- marker 支持；
- graph 连贯；
- NMF programme coherent；
- seed 稳定；
- 不被 stress/ribo/cell-cycle 主导。

## Program NMF

默认用 FastCNMF profile：`programs.method=fastcnmf`，`max_iter=50`，20 个 seed，K sweep 为 `[5, 8, 10, 12, 15, 20]`。`max_iter=50` 是当前内部 S2 和空间 S1 基准中更稳妥的默认值；如果 programme 不稳定，优先按 broad lineage 分层后重跑，而不是直接把 K 或迭代数大幅加高。

## Stress / ribosomal / proliferation

默认不删除这些细胞。先计算 sample 内 z-score/MAD outlier，再看是否主导 PC、是否跨 lineage 聚集、parent identity 是否清楚、是否与 doublet/ambient/low counts/high mt/bad FOV 共现。

决策顺序：

1. 保留 biology embedding；
2. 如影响 cell type clustering，从 identity HVG 排除相关 genes；
3. 或只对 identity embedding 做 regression；
4. 如果 parent identity 清楚，保留为 cell_state；
5. 只有多证据 artifact 才 remove。

## RCTD-py mode

低分辨率 spot/ROI/large bin 默认 `full`，因为它估计所有 cell types 的连续 mixture。`multi` 用于稀疏解释和敏感性分析，尤其是每个空间单位预期只含少数主导细胞类型时。

## CCC

FastCCC 是筛选，不是机制证明。复杂 receptor、integrin、ECM、MHC、cytokine receptor 等必须经过 complex-aware validation。空间项目还要检查 sender 和 receiver 是否相邻或共定位。
