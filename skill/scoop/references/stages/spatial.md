# Spatial views, graphs and RCTD-py


Inputs:

- Spatial H5AD view.
- Annotated single-cell reference H5AD.
- Spatial coordinates and optional SpatialData sidecar.

Script:

```bash
scoop spatial --config runs/<run_id>/config/run.yaml --dry-run
scoop spatial --config runs/<run_id>/config/run.yaml
```

Default rule:

- Low-resolution spot/ROI/large bin: RCTD-py `full` mode.
- Sparse or near-single-cell units: run `multi` or `doublet` only when justified.
- Check the specific backend with a representative real pilot, not only a
  `spatial: ok` flag. Use a recorded project adapter if the CLI/backend is unusable.

Sidecar outputs:

- `deconvolution/rctd_command.json`.
- RCTD-py output H5AD or tables, depending on configured command.
- Deconvolution weights should be externalised as Parquet.

Review triggers:

- `full` and `multi` disagree strongly.
- Weights are diffuse.
- Reference lacks expected tissue cell types.
- Spatial domains align with bad FOV/tile artefacts.

Before clinical interpretation check reference coverage/distinguishability, unit
resolution, patient/outcome links and tissue sampling. Raw weights and normalized
proportions stay separate; raw sums need not be one. Convergence/parity do not
prove biology. Bins/ROIs are not patients; deconvolution does not prove contact,
absolute counts or function. Review biology separately from numerics; unresolved
clinical mapping need not block descriptive work.
