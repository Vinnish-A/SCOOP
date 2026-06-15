# Spatial views, graphs and RCTD-py


Inputs:

- Spatial H5AD view.
- Annotated single-cell reference H5AD.
- Spatial coordinates and optional SpatialData sidecar.

Script:

```bash
python scripts/05_spatial_rctd.py --config runs/<run_id>/config/run.yaml --dry-run
python scripts/05_spatial_rctd.py --config runs/<run_id>/config/run.yaml
```

Default rule:

- Low-resolution spot/ROI/large bin: RCTD-py `full` mode.
- Sparse or near-single-cell units: run `multi` or `doublet` only when justified.

Sidecar outputs:

- `deconvolution/rctd_command.json`.
- RCTD-py output H5AD or tables, depending on configured command.
- Deconvolution weights should be externalised as Parquet.

Review triggers:

- `full` and `multi` disagree strongly.
- Weights are diffuse.
- Reference lacks expected tissue cell types.
- Spatial domains align with bad FOV/tile artefacts.
