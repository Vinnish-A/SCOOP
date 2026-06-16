# Scripts

这些脚本是 Agent 的执行接口。每个脚本只执行一个模块，并且通过 `--config` 读取 run 配置。

执行顺序：

```bash
python scripts/00_validate_input.py --config runs/<run_id>/config/run.yaml
python scripts/01_qc_scrublet.py --config runs/<run_id>/config/run.yaml
python scripts/02_core_analysis.py --config runs/<run_id>/config/run.yaml
python scripts/03_fast_consensus_nmf.py --config runs/<run_id>/config/run.yaml
python scripts/04_annotation_markers.py --config runs/<run_id>/config/run.yaml
python scripts/04b_tumor_fastcnvpy.py --config runs/<run_id>/config/run.yaml --gene-metadata data/external/references/gene_metadata.tsv
python scripts/04c_annotation_decide.py --config runs/<run_id>/config/run.yaml
python scripts/04d_annotation_commit.py --config runs/<run_id>/config/run.yaml --decisions runs/<run_id>/04_annotation/decisions/annotation_decision_template.json
python scripts/05_spatial_rctd.py --config runs/<run_id>/config/run.yaml --dry-run
python scripts/06_ccc_fastccc.py --config runs/<run_id>/config/run.yaml --dry-run
python scripts/08_prepare_pseudobulk.py --config runs/<run_id>/config/run.yaml
fastde deseq2 runs/<run_id>/07_de/pseudobulk/<cell_type> condition ctrl test
python scripts/07_prune_and_manifest.py --config runs/<run_id>/config/run.yaml
```

空间和 CCC 模块默认在配置里关闭。需要使用时，把 `spatial.enabled` 或 `ccc.enabled` 改成 `true`，并补齐 reference、LR resource、CellPhoneDB database 等路径。

`02_core_analysis.py` 现在调用 FastCore runner。SOP 入口仍是 `02_core`，但计算内核由 `core.engine` 控制；默认 `fastcore` 会先做 capability planning，FastCore 后端不可用或未通过策略选择时整体回退到唯一的 `scanpy_legacy` 后端。选中 `fastcore_oom` 时，脚本不会预读 H5AD，而是把路径直接交给 OOM backend。

FastCore backend 基准入口：

```bash
python scripts/fastcore/benchmark_core_backends.py \
  --config runs/<run_id>/config/run.yaml \
  --input runs/<run_id>/artifacts/adata_qc.h5ad
```

`--dry-run` 会输出外部命令而不执行，适合 Agent 在正式运行前检查参数、路径和 mode 选择。

DE 分成两类，不要混用：

- marker genes：回答“谁定义这个细胞群”，默认用 `fastde markers --method cosg`，也保留 Wilcoxon 备选；
- condition DE：回答“同一细胞群在条件改变后发生了什么”，默认用 `fastde deseq2` 的 pseudobulk DESeq2-like NB Wald。
- abundance：回答“样本级细胞类型/状态比例是否随条件、结局或生存改变”，默认用 `fastde abundance`，不做 per-cell p value。

需要和 R DESeq2 reference 对照时：

```bash
Rscript r/run_pseudobulk_deseq2.R runs/<run_id>/07_de/pseudobulk/<cell_type> condition ctrl test
python scripts/fastde/benchmark_deseq2.py --output-dir tmp/fastde_deseq2_benchmark
```

```bash
python scripts/fastde/benchmark_abundance.py --output-dir tmp/fastde_abundance_benchmark
```

legacy logCPM/Welch 和 edgeR QL 对照仍保留，主要用于回归检查：

```bash
Rscript r/run_pseudobulk_edger.R runs/<run_id>/07_de/pseudobulk/<cell_type> condition ctrl test
python scripts/benchmark_pseudobulk_de_python_vs_edger.py --output-dir tmp/pseudobulk_de_python_benchmark
```
