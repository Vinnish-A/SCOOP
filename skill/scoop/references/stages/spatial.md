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

The managed default now calls the pinned RCTD-py 0.3.9 API directly; the old
incompatible command template is removed from new configs. Explicit external
adapters remain possible with `spatial.rctd.adapter: external` and a supplied
`command_template`. `spatial.unit_type` controls the mode decision.

Set `spatial.rctd.q_matrices` or `RCTD_Q_MATRICES` to a staged upstream lookup
file; otherwise the upstream package uses its documented cache/download path.
Use native full reference counts before gene intersection. `exclude_genes`
removes held-out genes from reference fitting while preserving full spatial UMI
exposure. Save raw weights separately from nonnegative normalized visualizations.
Multi mode reconstructs the sparse result from `sub_weights/cell_type_indices`;
`weights` is its initial full fit. Doublet mode similarly keeps pair weights and
spot-class diagnostics. These estimates do not commit biological cell labels.
