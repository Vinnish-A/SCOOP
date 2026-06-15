# 01. 项目结构

SCOOP（Single Cell Omics Operating Protocol）的项目结构按“原始数据、运行结果、最终交付”分离。

```text
SCOOP/
  project.yaml
  environment.yml
  environment_omicverse.yml
  configs/
    default_run.yaml
    h5ad_schema.yaml
    omicverse_reuse_policy.yaml
  data/
    raw/
    external/
      references/
      lr_resources/
    processed/
  h5ad/
    canonical/
    spatial_views/
  spatialdata/
  runs/
    template/
      config/run.yaml
      01_qc/
      02_core/
      03_programs/
      04_annotation/
      05_spatial/
      06_ccc/
      07_de/
      reports/
      artifacts/
  src/scsp_agent_sop/
  scripts/
  docs/
```

## `data/`

`data/raw/` 是只读区。原始矩阵、图像、transcript coordinates 和 metadata 不在这里被修改。任何清洗或转换写入 `data/processed/` 或 run 目录。

`data/external/` 保存 reference atlas、marker database、ligand-receptor resource、CellPhoneDB database、gene set 文件等。外部资源需要版本和 hash。

## `h5ad/`

`h5ad/canonical/` 保存主对象，例如输入、QC 后、注释后和最终 H5AD。`h5ad/spatial_views/` 保存 spot、bin、cell segmentation 等不同空间分辨率的 H5AD view。

## `runs/`

每次分析一个 run。run 内部的七个模块目录互相隔离，大型结果表不写入 H5AD。run 完成后，`artifacts/final_adata.h5ad` 是瘦身后的交付对象，`reports/` 是人类阅读结果，`logs/decision_log.jsonl` 是审计轨迹。

## 为什么不把所有东西放进一个结果文件

FastCCC、CellPhoneDB、DE、marker、NMF weights、Leiden sweep 都是长表。把这些塞进 `adata.uns` 会导致 H5AD 难以加载、难以 diff、难以审计，也很容易破坏后续分析。外置表格更适合查询、版本控制和统计复核。
