# 04. Agent Execution Interfaces

Each script accepts `--config` and reads `runs/<run_id>/config/run.yaml` by default. Each script does one thing and, on completion:

1. Writes result tables;
2. Updates lightweight H5AD state;
3. Updates `file_registry`;
4. Appends to `logs/decision_log.jsonl`.

## QC

```bash
scoop qc --config runs/<run_id>/config/run.yaml
```

Key parameters are in the configuration file’s `qc` section:

- `counts_layer`
- `threshold_method`
- `mad.low_counts`
- `mad.high_mt`
- `scrublet.expected_doublet_rate`

Tuning rule: if more than 30% of a sample fails QC, do not simply widen the thresholds. First inspect the joint distributions of counts, genes, mt, ribo, hb, ambient, and doublets in that sample, then decide whether to enable GMM fallback or manual review.

## Core analysis

```bash
scoop core --config runs/<run_id>/config/run.yaml
```

Key parameters:

- `core.n_top_hvg`
- `core.n_pcs`
- `core.neighbors_n_neighbors`
- `core.leiden_resolutions`
- `core.batch_correction.method`

Tuning rule: initially retain `n_neighbors=15` and `n_pcs=30` by default. Adjust these parameters only when rare populations disappear, continuous trajectory structure breaks, or markers are lost after batch correction.

## FastCNMF programme discovery

```bash
scoop programs --config runs/<run_id>/config/run.yaml
```

Key parameters:

- `programs.method`
- `programs.k_grid`
- `programs.seeds`
- `programs.max_iter`
- `programs.stability_threshold`
- `programs.n_top_hvg`

Defaults are `programs.method=fastcnmf`, `programs.max_iter=50`, and 20 replicate seeds. Tuning rule: if fewer stable programmes than expected are found, first restrict the run to a broad lineage. Do not simply increase K; a larger K may only split ribosomal/stress/cell-cycle programmes.

## Spatial RCTD-py

```bash
scoop spatial --config runs/<run_id>/config/run.yaml
```

Key parameters:

- `spatial.rctd.mode_low_resolution`
- `spatial.rctd.mode_sparse_sensitivity`
- `spatial.rctd.expected_cells_per_unit`
- `spatial.rctd.command_template`

Use `full` for low resolution by default. When each spot/bin is expected to contain only a few cell types, or `full` weights are too diffuse, also run `multi` for sensitivity analysis.

## CCC

```bash
scoop ccc --config runs/<run_id>/config/run.yaml
```

Key parameters:

- `ccc.groupby`
- `ccc.lr_resource`
- `ccc.fastccc_command_template`
- `ccc.complex_validation.enabled`

Tuning rule: if `cell_type_lvl3` is too fine and leaves too few cells in some groups, fall back to `cell_type_lvl2`. Multimeric LRIs or LRIs central to a mechanism must pass the CellPhoneDB / LIANA gate.

## Prune and manifest

```bash
scoop deliver --config runs/<run_id>/config/run.yaml
```

This step removes temporary columns according to `configs/h5ad_schema.yaml`, writes `manifest.json`, and generates a slimmed-down H5AD.
