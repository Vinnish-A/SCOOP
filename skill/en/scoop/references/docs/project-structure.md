# 01. Project Structure

As a software repository, SCOOP (Single Cell Omics Operating Protocol) stores only code, configuration templates, documentation, tests, and run templates. Real project data and actual run outputs belong to the calling project or local test workspace and do not enter the SCOOP git repository.

## SCOOP repository

```text
SCOOP/
  configs/
    default_run.yaml
    h5ad_schema.yaml
    omicverse_reuse_policy.yaml
  templates/
    run/
      config/run.yaml
      qc/
      core/
      programs/
      annotation/
      spatial/
      ccc/
      differential-expression/
      reports/
      artifacts/
  src/scsp_agent_sop/
  src/scoop/
  scripts/
    r/
  docs/
  tests/
```

`scripts/` stores command-line entry points that agents and analysts can execute directly. Python SOP scripts live at the root of `scripts/` or in the corresponding fast subdirectory; R reference/regression scripts live in `scripts/r/` rather than separately occupying the project root.

`templates/run/` is the run-directory template for new projects. To use it, copy it to `runs/<run_id>/` in the calling project or to the local `.scoop_local/runs/<run_id>/`.

## Calling project

When SCOOP serves other projects as an MCP service or command-line tool, the calling project owns its data and run directories:

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

`data/raw/` is read-only. Original matrices, images, transcript coordinates, and metadata are not modified here. Write any cleaning or transformations to `data/processed/` or the run directory.

`data/external/` stores reference atlases, marker databases, ligand-receptor resources, CellPhoneDB databases, gene-set files, and similar resources. External resources require versions and hashes.

`data/h5ad/` stores project-level H5AD objects and quick-test H5ADs. `data/h5ad/canonical/` is the main object area, for example for input, post-QC, post-annotation, and final H5ADs; `data/h5ad/spatial_views/` stores H5AD views at different spatial resolutions, such as spots, bins, and cell segmentations.

`data/spatialdata/` stores SpatialData/Zarr, GeoParquet, images, polygons, transcript coordinates, and spatial coordinate-transformation sidecars. It is separate from H5AD to avoid stuffing large image and geometry objects into H5AD.

## Local test workspace

Local tests in the SCOOP repository may use `.scoop_local/`, which is ignored by `.gitignore`:

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

Store actual quick-test H5ADs, private data, raw spatial samples, and benchmark outputs in `.scoop_local/` or `tmp/`; do not commit them to git.

## `runs/`

Use one run per analysis. Module directories within a run are isolated from each other, and large result tables are not written into H5AD. On completion, `artifacts/final_adata.h5ad` is the slimmed-down delivery object, `reports/` holds human-readable results, and `logs/decision_log.jsonl` is the audit trail.

## Why not put everything into a single result file

FastCCC, CellPhoneDB, DE, markers, NMF weights, and Leiden sweeps are all long tables. Stuffing them into `adata.uns` makes H5AD difficult to load, diff, and audit, and can easily break downstream analysis. External tables are better suited to querying, version control, and statistical review.
