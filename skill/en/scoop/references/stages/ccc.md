# FastCCC and complex-aware validation


Inputs:

- Annotated H5AD.
- Complex-aware ligand-receptor resource.
- Optional CellPhoneDB database zip.

Script:

```bash
scoop ccc --config runs/<run_id>/config/run.yaml --dry-run
scoop ccc --config runs/<run_id>/config/run.yaml
```

Optional validation:

```bash
scoop ccc --config runs/<run_id>/config/run.yaml --validate-cellphonedb
scoop ccc --config runs/<run_id>/config/run.yaml --validate-liana
```

Default tool:

- FastCCC for primary screening. The default input is `artifacts/adata_annotation_committed.h5ad`; only cells with
  `annotation_status` `accepted`/`labeled` in groups of at least `ccc.min_cells_per_group` are tested. When cells
  were dropped, FastCCC reads the filtered `artifacts/adata_ccc_indexed.h5ad`, and the receipt lists both.

OmicVerse reuse:

- `omicverse.single.run_cellphonedb_v5` for complex-sensitive interactions.
- `omicverse.single.run_liana` for rank-aggregate validation.

Sidecar outputs:

- `fastccc/tables/complex_sensitive_lris.parquet`.
- FastCCC output directory from configured command.
- CellPhoneDB/LIANA validation tables.

Review triggers:

- Any mechanistic claim depends on one LRI.
- LRI contains multimeric ligand or receptor.
- Spatial claim lacks adjacency support.
