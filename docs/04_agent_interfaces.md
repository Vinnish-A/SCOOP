# 04. Agent 执行接口

每个脚本都接受 `--config`，默认读取 `runs/<run_id>/config/run.yaml`。脚本只做一件事，并且在完成时：

1. 写出结果表；
2. 更新 H5AD 轻量状态；
3. 更新 `file_registry`；
4. 追加 `logs/decision_log.jsonl`。

## QC

```bash
python scripts/01_qc_scrublet.py --config runs/<run_id>/config/run.yaml
```

关键参数在配置文件的 `qc` 段：

- `counts_layer`
- `threshold_method`
- `mad.low_counts`
- `mad.high_mt`
- `scrublet.expected_doublet_rate`

调参规则：如果某样本 QC fail 比例大于 30%，不要直接加宽阈值。先检查该样本的 counts、genes、mt、ribo、hb、ambient、doublet 共同分布，再决定是否启用 GMM fallback 或人工 review。

## Core analysis

```bash
python scripts/02_core_analysis.py --config runs/<run_id>/config/run.yaml
```

关键参数：

- `core.n_top_hvg`
- `core.n_pcs`
- `core.neighbors_n_neighbors`
- `core.leiden_resolutions`
- `core.batch_correction.method`

调参规则：默认先保持 `n_neighbors=15` 和 `n_pcs=30`。只有在 rare populations 消失、trajectory 连续结构断裂、或 batch correction 后 marker 丢失时，才调整这些参数。

## FastCNMF programme discovery

```bash
python scripts/03_fast_consensus_nmf.py --config runs/<run_id>/config/run.yaml
```

关键参数：

- `programs.method`
- `programs.k_grid`
- `programs.seeds`
- `programs.max_iter`
- `programs.stability_threshold`
- `programs.n_top_hvg`

默认 `programs.method=fastcnmf`，`programs.max_iter=50`，20 个 replicate seeds。调参规则：如果稳定 programme 少于预期，先缩小到 broad lineage 内运行。不要简单扩大 K；K 增大可能只会拆分 ribosomal/stress/cell-cycle 程序。

## Spatial RCTD-py

```bash
python scripts/04_spatial_rctd.py --config runs/<run_id>/config/run.yaml
```

关键参数：

- `spatial.rctd.mode_low_resolution`
- `spatial.rctd.mode_sparse_sensitivity`
- `spatial.rctd.expected_cells_per_unit`
- `spatial.rctd.command_template`

默认低分辨率用 `full`。当每个 spot/bin 预期只有少数 cell types，或者 `full` 权重过于 diffuse，再运行 `multi` 做敏感性分析。

## CCC

```bash
python scripts/05_ccc_fastccc.py --config runs/<run_id>/config/run.yaml
```

关键参数：

- `ccc.groupby`
- `ccc.lr_resource`
- `ccc.fastccc_command_template`
- `ccc.complex_validation.enabled`

调参规则：如果 `cell_type_lvl3` 太细导致某些 group 细胞太少，回退到 `cell_type_lvl2`。多聚体 LRI 或机制关键 LRI 必须过 CellPhoneDB / LIANA 守门。

## Prune and manifest

```bash
python scripts/06_prune_and_manifest.py --config runs/<run_id>/config/run.yaml
```

该步骤会基于 `configs/h5ad_schema.yaml` 删除临时列、写出 `manifest.json`、生成瘦身 H5AD。
