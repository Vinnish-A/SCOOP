# FastCCC and complex-aware validation


Inputs:

- Annotated H5AD.
- Complex-aware ligand-receptor resource.
- Optional CellPhoneDB database zip.

Script:

```bash
python scripts/06_ccc_fastccc.py --config runs/<run_id>/config/run.yaml --dry-run
python scripts/06_ccc_fastccc.py --config runs/<run_id>/config/run.yaml
```

Optional validation:

```bash
python scripts/06_ccc_fastccc.py --config runs/<run_id>/config/run.yaml --validate-cellphonedb
python scripts/06_ccc_fastccc.py --config runs/<run_id>/config/run.yaml --validate-liana
```

Default tool:

- FastCCC for primary screening.

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
