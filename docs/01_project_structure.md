# 01. 项目结构

SCOOP（Single Cell Omics Operating Protocol）作为软件仓库时只保存代码、配置模板、文档、测试和 run 模板。真实项目数据和真实运行输出属于调用方项目或本地测试工作区，不进入 SCOOP git 仓库。

## SCOOP 仓库

```text
SCOOP/
  configs/
    default_run.yaml
    h5ad_schema.yaml
    omicverse_reuse_policy.yaml
  templates/
    run/
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
  src/scoop/
  scripts/
    r/
  docs/
  tests/
```

`scripts/` 保存 Agent 和分析人员可直接执行的命令行入口。Python SOP 脚本放在 `scripts/` 根或对应 fast 子目录；R 参考/回归脚本放在 `scripts/r/`，不再单独占用项目根目录。

`templates/run/` 是新项目 run 目录模板。使用时复制到调用方项目的 `runs/<run_id>/`，或复制到本地 `.scoop_local/runs/<run_id>/`。

## 调用方项目

SCOOP 作为 MCP 服务或命令行工具服务其它项目时，调用方项目拥有自己的数据和运行目录：

```text
Project/
  project.yaml
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
    <run_id>/
      config/run.yaml
      artifacts/
      reports/
      logs/
```

`data/raw/` 是只读区。原始矩阵、图像、transcript coordinates 和 metadata 不在这里被修改。任何清洗或转换写入 `data/processed/` 或 run 目录。

`data/external/` 保存 reference atlas、marker database、ligand-receptor resource、CellPhoneDB database、gene set 文件等。外部资源需要版本和 hash。

`data/h5ad/` 保存项目级 H5AD 对象和快速测试 H5AD。`data/h5ad/canonical/` 是主对象区，例如输入、QC 后、注释后和最终 H5AD；`data/h5ad/spatial_views/` 保存 spot、bin、cell segmentation 等不同空间分辨率的 H5AD view。

`data/spatialdata/` 保存 SpatialData/Zarr、GeoParquet、图像、多边形、transcript coordinates 和空间坐标转换 sidecar。它和 H5AD 分开，是为了避免把图像和几何大对象塞进 H5AD。

## 本地测试工作区

SCOOP 仓库本地测试可以使用 `.scoop_local/`，该目录被 `.gitignore` 忽略：

```text
.scoop_local/
  data/
    h5ad/
    raw/
    external/
    spatialdata/
  runs/
  tmp/
```

真实快速测试 H5AD、私有数据、空间原始样本和 benchmark 输出都放在 `.scoop_local/` 或 `tmp/`，不提交到 git。

## `runs/`

每次分析一个 run。run 内部的模块目录互相隔离，大型结果表不写入 H5AD。run 完成后，`artifacts/final_adata.h5ad` 是瘦身后的交付对象，`reports/` 是人类阅读结果，`logs/decision_log.jsonl` 是审计轨迹。

## 为什么不把所有东西放进一个结果文件

FastCCC、CellPhoneDB、DE、marker、NMF weights、Leiden sweep 都是长表。把这些塞进 `adata.uns` 会导致 H5AD 难以加载、难以 diff、难以审计，也很容易破坏后续分析。外置表格更适合查询、版本控制和统计复核。
