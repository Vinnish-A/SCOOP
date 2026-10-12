# 实例二：BRCA / ARTEMIS 的相对比例、UCell 与 pCR

检索词：zzu_gp / relative fraction / parent lineage / UCell / rare state /
pCR / OR / LAM / LPL / forest / 临床注释。
核查日期：2026-10-07。复用已有患者表和回归结果，没有新拟合。

## 要回答的问题

患者中一个状态“占父群多少”，与这个状态的表达程序“有多强”分别怎样关联
新辅助治疗响应？数量和功能分开测量，才能从 LAM 的生态背景自然引出 LPL。

## 实际队列与两种测量

项目根目录：`/mnt/nfs_data/xzh/BRCA`。
QC 后 417,698 细胞、101 患者；原 49 状态森林有 82 位结局可用患者，
45 pCR、37 RD。固定原状态及已审核身份，不通过改标签获得临床方向。

旧方法来自 `zzu_gp/R/0Fig1.5.R:692–780`，在本项目具体落实为：

| 测量 | 实际定义 | 不能混称的对象 |
| --- | --- | --- |
| 常见状态比例 | 每患者 state_count / 对应父谱系 count | 全样本总细胞比例 |
| 稀有状态 UCell | 固定发表 marker；全基因排名；每患者全部细胞中严格 >Q70 的分数均值 | 阳性细胞比例、counts pseudobulk、LAM 实际数量 |
| LAM 计数对照 | 原 C1Q/CD36 状态 ∩ approved Myeloid，除以同范围所有原 mac-* 细胞 | 全患者 UCell 程序强度 |
| 髓系 LPL 输出 | approved Myeloid 内 counts pseudobulk 的 log2(1+CPM) | LAM 专属表达、UCell、细胞比例 |

旧稀有规则是“检出该状态患者的全细胞比例均值 <0.015”，不是全队列零补齐均值。
C1Q=0.013348、CD36=0.007004，故旧森林两行都走 UCell。
这属于项目方法记录，不要求所有新项目沿用该门槛。

实际 UCell 在完整 32,354 基因排名，maxRank=1500、average ties、missing genes
impute；49 个固定签名含发表的 14–20 marker，两 LAM 各 20。来源是 canonical
`yan_counts.h5ad` 的 counts 层，不是仅 HVG 的整合矩阵。主汇总在患者全部细胞，
父谱系内上尾另存。已有 R UCell 2.14.0 对 72 个实际细胞 ×49 签名的验证，最大
绝对误差 2.220446e−16。这是历史实现等价核对，不意味着不同 marker/scope 等价。

常见状态按父群有无决定覆盖：父群存在而目标缺失填 0，父群不存在比例为 NA。
因此旧森林常见状态的患者数有 T=82、Mye=79、B=70、Fibro=20，不能统一写 82。

## 已有临床结果

森林统一 pCR=1/RD=0，连续暴露按模型患者内 1 SD。
原输出若以 RD 为事件，取倒数得到 OR_pCR 并互换 CI 边界，P 不变。

| 实际指标 | 患者数 | 模型 | OR_pCR | 95% CI | nominal P |
| --- | ---: | --- | ---: | --- | ---: |
| C1Q 全患者 UCell 上尾均值 | 82 | unadjusted | 0.7662 | 0.4921–1.1931 | 0.2385 |
| CD36 全患者 UCell 上尾均值 | 82 | unadjusted | 0.8587 | 0.5546–1.3296 | 0.4947 |
| C1Q / Mac 原始比例 | 70 | unadjusted | 1.7458 | 1.0415–2.9265 | 0.03449 |
| C1Q / Mac 原始比例 | 70 | composition | 1.6866 | 0.9774–2.9101 | 0.06038 |
| 髓系 LPL log2(1+CPM) | 70 | composition | 0.5025 | 0.2902–0.8702 | 0.01404 |

82 人 UCell 与 70 人比例不是同一个暴露/覆盖。后两类的同 70 人对照包括
37 pCR、33 RD；定义父群的全集共 22,744 个合格巨噬细胞，逐患者要求父群 ≥30。
稀有 LAM 分子本身没有 ≥30 门槛；CD36 仅 8 位非零、62 位真零。
composition 协变量固定为 Myeloid、T、Epithelial 三个全细胞比例。

实际可讲的结论是：C1Q 捕获比例与髓系 LPL 输出的临床方向不同；本数据中
LPL 输出的响应关联更清楚。UCell 没有得到 LAM 总体程序强度的显著响应关联，
不能换称为已显著的 LAM 比例。固定 logit 敏感性是另一条测量，不能拿其 P
替代原始比例的 P；相同患者的另一模型也不是独立队列验证。

## 直接可用的表与来源

- [全部 49 状态主森林表](data/brca-state-primary-or.tsv)：带原状态编号、测量类型、
  OR/CI/P、覆盖、原 RD 方向和 q；不只保留 LAM 两行。
- [同 70 人数量/功能对照](data/brca-lam-lpl-comparison.tsv)：全部 10 行原对照。
- [LAM 稀有方法审查](data/brca-lam-rare-audit.tsv)：实际细胞数、比例、签名来源。
- 原 49 状态结果：`results/fig1_vstory_20261005/rebuild_scrna/`，其中
  `forest49_primary_OR_pCR.tsv`、`analysis_config.json`、
  `REBUILD_SCRNA_METHODS.zh.md` 保存主结果和参数。
- 原患者计数表：`results/LAM_outcome_EDA_20261005/SC/patient_fixed_LAM_allmac_estimands.tsv`；
  原模型表与报告同目录，代码位于 `analysis/LAM_outcome_EDA_20261005/SC/`。

数值副本与源路径/哈希见 [案例清单](data/source-manifest.json)。

## 如何做成信息密集的主图

以状态身份/功能热图与临床森林为共同索引，常见比例和稀有 score 语义分组，
加入窄患者轨和临床指标对照。在森林前以脂质程序成员和谱系表达引出 LPL，
侧边用同患者数量–功能比较解释为何聚焦 LPL；不要仅为显著性跳换状态名称。
单细胞结局接原生组织、生态位和患者内空间对照，使用真实配对 crosswalk。
热图和森林实际数据边缘对齐，临床轨共享患者顺序，图例插入空隙。
实验槽留给 LAM/LPL 定位和功能验证，不用计算 score 假充实验读数。

新项目从冻结身份和完整基因表重建 patient×state 表，再按实际 endpoint
输出 OR 或 HR。旧 0.015/Q70/top20 参数是可审阅的范例；改变父群、marker 或
汇总后，给新测量独立名称和完整结果，不声称只是换了同一个模型的计算后端。
算法参照 [UCell 官方仓库](https://github.com/carmonalab/UCell)。
