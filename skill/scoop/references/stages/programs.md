# NMF programme discovery

For a formal fine-subtype analysis, first accept a biologically meaningful broad
parent class through review of clustering, conventional DEG and tissue markers;
numeric clusters or unreviewed merged ranges do not by themselves define it.
Broad-class acceptance does not require NMF. After that acceptance, parent NMF
and preparatory within-parent clustering/DEG may run in parallel. Formal subtype
analysis requires the parent and technical prerequisites to be complete, including
the appropriately batch-handled parent NMF; formal subtype adjudication requires
all evidence in the [annotation workflow](../markers/annotation-workflow.md).
Retain the declared UID/source scope and mixed or competing contexts. An exploratory
programme on an unreviewed range is not a completed biological parent analysis.
Programme argmax is not a cell identity classifier. If parent evidence is missing,
record evidence-pending progress and actively obtain it; only adequate but
conflicting/mixed reviewed evidence supports an unresolved identity.

For cross-dataset discovery, follow the
[metadata and batch contract](../docs/design-principles.md#元信息与批次契约).
Run and verify the supported expression-level or model-internal batch interface
before NMF, recording actual dataset keys, layer/scale and input/output bindings.
A latent-only Harmony result cannot satisfy this expression-processing step;
if the chosen route lacks the required interface, record that gap and use a
validated compatible route before calling the main analysis batch-handled.

Native consensus uses one-to-one component matches across independent random
seeds. Report per-component cross-seed cosine, the weakest component and actual
convergence; same-seed components and singleton groups cannot establish stability.
One seed means unassessed stability. A high aggregate median can hide a weak
programme and must not qualify every component. The native K rule requires the
minimum component score to pass; a fallback selection remains exploratory.

`programs.hvg_key` defaults to `highly_variable_biology`; set it explicitly to
null for an all-input-feature sensitivity. Native MT/HB exclusions still apply.
Record input/discovery feature counts and the discovery-axis hash separately
from any full-gene projection or association. The component evidence table and
selected-K convergence/refit diagnostics are part of the programme result.


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

For standalone `fast-prepare -> fast-factorize -> fast-consensus`, record the
actual nonnegative input scale, feature selection, normalization and any technical
correction; the argument name `corrected-h5ad` does not prove correction occurred.
The full official cNMF consensus path reads the normalized-count H5AD: prepare
with `--norm-store h5ad` or `both`. NPY-only preparation is appropriate only for
paths that actually support it; factorization success alone does not validate
full consensus. Preserve valid replicate files when repairing an output format,
verify the same normalized input/parameters and resume consensus instead of
refitting. Report exact/minibatch/GPU factorization backends separately.

Validate actual output axes and numeric semantics: consensus spectra/usage and
nonnegative expression spectra differ from signed gene association scores. Save
per-K independent seeds, iterations/convergence, stability, fixed-H usage and
sample/QC/rare-group coverage. A representative source-level probe establishes
backend execution, not programme completion for a newly integrated compartment.

For an anchor-fit route, freeze the complete eligible UID scope, then audit
parent-group × true sampling-unit coverage **after** assay/observability exclusions.
A balanced pool before exclusions can silently lose an eligible unit. Reuse valid
expression checkpoints, supplement only missing eligible units from their frozen
native caches, and preserve the original analysis. Disclose anchor versus full
scope counts; anchor fitting is not an all-cell discovery fit. Use native-depth
linear CP10k when that is the declared scale; do not renormalize the intersected
genes or call CP10k TPM.

Freeze anchor feature selection and training scales. Before projecting onto
unseen eligible rows, reproduce anchor usage with the actual official final
basis, normalization and training standard-deviation convention. Median
probability spectra, expression spectra and signed scores are not interchangeable
projection bases. Save this numeric proof and block-level UID/feature/parameter
checkpoints; unassessed assay rows remain unassessed rather than receiving
fabricated zero programme scores. Review per-K convergence and source/depth/QC
associations before biological naming or choosing a final K.

Report discovery HVGs separately from genes in the final expression/association
spectra. A full-gene refit output does not mean every gene participated in factor
discovery. Inspect whether identity and competing-lineage anchors were available
to discovery before interpreting a missing programme as biological absence.
High component silhouette and successful anchor reproduction establish numeric
properties, not biological specificity; neither justifies choosing K or assigning
cell identity by programme argmax. Preserve a weak biological result rather than
tuning K merely to reproduce the desired names.
