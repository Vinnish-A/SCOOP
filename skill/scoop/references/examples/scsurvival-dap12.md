# 实例一：DAP12 / GSE274546 的 scSurvival

检索词：scSurvival / OS / Cox / MIL / attention / cellular risk / 髓系 / MYD88。
核查日期：2026-10-07。这里读取已有真实输出，没有重新训练。

## 要回答的问题

髓系细胞内部的哪些表达状态承载患者生存风险？再看这些状态的身份、
MYD88 表达和患者捕获比例，判断“风险状态”与“某一群细胞变多”是否同向。
这是表达特征驱动的 scSurvival，不能以比例回归冒充。

## 实际输入与执行

项目根目录：`/mnt/sdb/xzh/Vproject/DAP12`。

| 对象 | 项目内实际路径 | 内容 |
| --- | --- | --- |
| 细胞特征 | `result/GSE274546_myeloid_survival/scsurvival_input/gse274546_myeloid_cell_features.csv.gz` | barcode、sample、patient、髓系状态、MYD88、Harmony 特征 |
| 临床表 | `result/GSE274546_myeloid_survival/os_primary_paired_samples.csv` | 56 个样本/患者；OStime、OS；46 个事件 |
| 原执行脚本 | `script/run_gse274546_scsurvival_myd88.py` | 直接调用原包 `scSurvivalRun` |
| 细胞输出 | `result/GSE274546_myeloid_survival/scsurvival/cell_scsurvival_scores_primary_os.csv` | 39,842 细胞，54 个实际覆盖患者 |
| 患者输出 | `result/GSE274546_myeloid_survival/scsurvival/sample_scsurvival_hazards_primary_os.csv` | 54 行、45 个事件；time/status/patient_hazards |
| 状态输出 | `result/GSE274546_myeloid_survival/scsurvival/cluster_scsurvival_myd88_summary.csv` | 15 个 cluster 的风险、attention 和表达汇总 |

实际输出缺少临床表中的 `P28T1`、`P75T1`，不能把 56 当模型输出覆盖人数。
本次只核对覆盖差异，没有推定缺失原因。

旧脚本选前 30 个 Harmony 维度，StandardScaler 后置于 `obsm['X_harmony']`；
使用 `feature_flavor='Custom'`、`feature_key='X_harmony'`。AnnData 的一列全零
X 是接口占位，实际输入是表达衍生的 Harmony 特征，不能说模型训练在零表达上。
参数包括 5 heads、dropout 0.3、epochs 240、patience 60、lr 0.001、
instance_batch_size 2000、sample_balance=True、predict_nMC=32、validate=False。
这些记录用于解释旧产物，不作为其他队列的固定最优参数。

## 已有数值长什么样

可直接读取 [完整 15 cluster 表](data/dap12-scsurvival-summary.csv)。
下面是其中三行；风险是旧表的 pooled-cell mean，不是临床 HR。

| 状态 | 细胞数 | 检出样本数 | mean hazard_adj | mean attention |
| --- | ---: | ---: | ---: | ---: |
| TAM_SELENOP_MRC1 | 2,779 | 45 | −2.5290 | 0.6277 |
| Myeloid_Inflammatory_CCL3_IL1A（cluster 5） | 2,522 | 51 | −2.3242 | 0.7479 |
| TAM_APOC2_MT3 | 150 | 30 | −5.0002 | 0.7982 |

旧状态表另接了一套比例模型的结果：SELENOP/MRC1 行的 OS HR=0.7435、
P=0.0430；CCL3/IL1A 行的 OS HR=0.7782、P=0.2620，但其
OR_no_later_recur=0.4614、P=0.0290。这些来自
`result/GSE274546_myeloid_survival/myeloid_scSurvival_like_results.csv`，不是
scSurvival 模型输出，也不能仅凭合表就声称使用了相同患者和效应单位。
下一次比较须从各模型输入/receipt 核对覆盖和单位。

## 可迁移的分析与图形动作

以细胞风险连续场为入口，接“患者 × 状态”的风险分布、attention、身份/功能
成员谱，最后在同一患者集比较状态比例与临床效应。风险分布、表达身份和
捕获数量各自承担一个问题；用共享状态编号和顺序建立视觉中心。
可将窄临床轨和比例效应贴在状态矩阵侧面，原组织图提供空间定位，实验槽留给
分选/定位与功能扰动。不要把 attention 再画成另一个未经定义的“细胞比例”。

这个例子留下三条直接经验：

- 低 OS 风险与少复发并不必然同向，终点必须分别读。
- APOC2_MT3 的风险很低，但 top markers 以线粒体基因为主；需结合质量和身份
  解释，不能靠最低风险自动命名为“保护性巨噬细胞”。
- `validate=False` 的旧输出是拟合风险画像。新项目若声称预测能力，应增加
  患者级分割的 held-out 风险和临床基线；不把旧拟合值当外部验证。

旧脚本还有风险分位数、比例方向及 IFN 例外的选群规则；它是历史筛选记录，
不应成为 Scoop 的自动生物学注释规则。复用细胞风险→状态功能→患者对照的
组织方式，重新确定本项目的状态解释和评价。

来源与哈希见 [案例清单](data/source-manifest.json)。算法依据为
[原 scSurvival 仓库](https://github.com/cliffren/scSurvival)。
