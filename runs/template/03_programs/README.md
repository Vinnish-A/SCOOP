# FastCNMF programme discovery


Inputs:

- `artifacts/adata_core.h5ad`.
- `layers["log1p_norm"]`.
- `var["highly_variable_biology"]`.

Script:

```bash
python scripts/03_fast_consensus_nmf.py --config runs/<run_id>/config/run.yaml
```

Default tool:

- FastCNMF profile with exact coordinate-descent NMF, K sweep, 20 replicate seeds and `max_iter=50`.
- Existing configs that use `fast_consensus_nmf` or `sklearn_nmf` remain supported as compatibility method names.

Optional validation:

```bash
python scripts/03_fast_consensus_nmf.py --config runs/<run_id>/config/run.yaml --validate-with-omicverse-cnmf
```

This calls `omicverse.single.cNMF` when installed. It is a validation fallback,
not the preferred programme discovery method.

H5AD writes:

- `obsm["X_nmf_usage"]`.
- `obs["dominant_nmf_program"]`.
- `obs["nmf_program_entropy"]`.

Sidecar outputs:

- `tables/nmf_k_sweep.tsv`.
- `tables/nmf_gene_weights.parquet`.
- `tables/nmf_usage.parquet`.

Review triggers:

- Programme stability below 0.70.
- A programme supports a key biological claim.
- A programme is sample-specific or QC-associated.
