# Scripts

这些脚本是 Agent 的执行接口。每个脚本只执行一个模块，并且通过 `--config` 读取 run 配置。

执行顺序：

```bash
python scripts/00_validate_input.py --config runs/<run_id>/config/run.yaml
python scripts/01_qc_scrublet.py --config runs/<run_id>/config/run.yaml
python scripts/02_core_analysis.py --config runs/<run_id>/config/run.yaml
python scripts/03_fast_consensus_nmf.py --config runs/<run_id>/config/run.yaml
python scripts/04_annotation_markers.py --config runs/<run_id>/config/run.yaml
python scripts/05_spatial_rctd.py --config runs/<run_id>/config/run.yaml --dry-run
python scripts/06_ccc_fastccc.py --config runs/<run_id>/config/run.yaml --dry-run
python scripts/08_prepare_pseudobulk.py --config runs/<run_id>/config/run.yaml
python scripts/07_prune_and_manifest.py --config runs/<run_id>/config/run.yaml
```

空间和 CCC 模块默认在配置里关闭。需要使用时，把 `spatial.enabled` 或 `ccc.enabled` 改成 `true`，并补齐 reference、LR resource、CellPhoneDB database 等路径。

`--dry-run` 会输出外部命令而不执行，适合 Agent 在正式运行前检查参数、路径和 mode 选择。
