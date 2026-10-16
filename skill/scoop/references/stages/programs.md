# NMF programme discovery


Inputs:

- `artifacts/adata_core.h5ad`.
- `layers["log1p_norm"]`.
- `var["highly_variable_biology"]`.

Script:

```bash
scoop programs --config runs/<run_id>/config/run.yaml
```

Default tool:

- Native `sklearn_multistart_nmf`: coordinate-descent NMF, K sweep, 20 independent random starts, component consensus and fixed-component usage refit; `max_iter=50` is the initial iteration budget.
- Existing configs that use `fast_consensus_nmf` or `sklearn_nmf` remain supported as compatibility method names.

Optional validation:

```bash
scoop programs --config runs/<run_id>/config/run.yaml --validate-with-omicverse-cnmf
```

This calls `omicverse.single.cNMF` when installed. It is an explicit validation path,
not an automatic response to low stability.

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


Native `scoop programs` reports `sklearn_multistart_nmf` explicitly. Older method
names are compatibility aliases, not a claim that the standalone FastCNMF engine
ran. Replicates use independent random initializations with BLAS threads bounded
across concurrent jobs. Delivered spectra are component consensus estimates and
usage is refitted with fixed consensus components. Report per-K convergence
fractions, iterations, seeds and actual method. A 50-iteration budget remains a
fast initial pass, not proof of convergence. One seed cannot establish stability.

OmicVerse validation is opt-in (`--validate-with-omicverse-cnmf` or
`programs.validate_with_omicverse_cnmf: true`); low stability records review status
without silently starting a more expensive algorithm.

The standalone `fastcnmf` prepare/factorize/consensus commands share the managed
Python environment, including pinned official `cnmf`. Invoke them via
`$SCOOP_HOME/bin/fastcnmf`; they expose cNMF-compatible full-consensus algorithms.
Choose this route explicitly when needed rather than treating native sklearn NMF
as the same estimator. No project-only cNMF installation is needed.
