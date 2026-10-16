# Large-atlas integration and reviewed supervision

Use scVI for unsupervised representations, then review identities before any
scPoli supervision. QC follows [QC](../stages/qc.md); cell labels follow the
[annotation workflow](../markers/annotation-workflow.md). Integrated expression
is not a substitute for raw-count pseudobulk in condition DE.

## Inputs and split

Keep integer `layers["counts"]`, unique cell/gene IDs, sample/library, donor,
study, condition, assay and tissue. Record duplicate-symbol handling and the
measured gene universe. Preserve legacy labels as audit metadata only.
Apply `final_use`; do not reintroduce confirmed doublets for batch mixing.
A lineage subset inherits original-library QC rather than refitting Scrublet.

Freeze train, validation and query by complete donor/study before evaluating
supervised mapping. A random cell split leaks patient-specific structure. Record
rare labels that cannot be independently held out.

## scVI execution

```bash
scoop atlas integrate --config run.yaml --max-epochs 60 --batch-size 2048 \
  --accelerator gpu --devices 1 --kl-warmup-epochs 20
```

These are example resource settings, not universal defaults. Native options are
available through `atlas.scvi.train_kwargs` and `atlas.scvi.plan_kwargs`, including
datasplitter/trainer settings. CLI overrides take precedence. Without an explicit
KL schedule, epoch warmup is one third of the training budget, capped at 400;
explicit epoch warmup must fit the budget. Step-based warmup remains explicit.

Use technical batch/covariates; do not automatically correct disease, treatment,
response or target identity. If study and disease are confounded, retain that
limitation. Record versions, model, batch size, accelerator/devices, train/validation
loss and actual training options. Model existence is not convergence.

The runtime builds the scVI graph and runs the shared parallel Leiden sweep.
`atlas.graph.leiden_resolutions`, `leiden_seeds` and `leiden_search` control its
grid/seeds/resources. By default the configured single resolution is checked at
three seeds. Keep per-job checkpoints and sweep/stability tables; a single seed
cannot establish stability. Check latent neighbours, marker preservation, sample
coverage and QC axes; UMAP mixing alone does not validate integration.

## Identity review and report

Review broad lineages first, then reconstruct within-lineage structure as needed.
Use DEG, positive/exclusion markers, tissue-matched literature, sample/donor/study
replication, nearest graph neighbours and programme evidence. Do not use automated
reference labels to assign identities. General stress/IFN/cycling differences
usually belong to state, not a new identity. Tumor status needs its own evidence.

Keep `cell_type_lvl1/2/3`, `cell_state`, `scvi_curated_cell_type`,
`annotation_confidence`, `annotation_status`, source and version distinct.
Report machine-readable marker/decision tables and an HTML report covering:

- counts and both conditional proportions by donor/study/condition;
- graph-neighbour differences, marker support/conflicts and merge/split rationale;
- programme identity/function/technical interpretations and their sample coverage;
- QC/doublet risks, unresolved alternatives and proposed supervision indices.

A cluster number or high programme score is not a label. High-confidence labels
need coherent marker/exclusion evidence; ambiguous, low/medium-confidence,
unreviewed legacy and suspicious cells remain unlabeled. Unresolved identities
can be reviewed without forcing a type.

## scPoli approval and evaluation

scPoli prototype training amplifies supervision errors. Require explicit current
approval of the annotation report and labels before training. Use `scoop atlas
approve` to bind the actual input H5AD, report, labels and exact labeled-index table
by SHA-256; `scoop atlas scpoli` rejects missing, stale or inconsistent gates.
Tool bypass does not replace this approval. Changing labels/scope invalidates it.

Only high-confidence reviewed cells supervise the model. If a coarse identity is
certain but a fine subtype is not, make a separate reviewed coarse decision;
do not promote uncertain fine labels. Use separate granularities when necessary.
Record supervised label sizes, donor coverage, prototype distances and uncertainty.

Evaluate complete held-out donors/studies, marker coherence, confusion, unknown
rate and drift. Calibrate distance thresholds on validation donors, then freeze
for queries; no universal prototype cutoff is implied. Reproducing training
labels is not independent validation. Report fitted and held-out results separately.

## Delivery

Keep counts, latent coordinates, graph descriptors, reviewed labels/confidence,
split and supervision flags. Keep models, training history, full evidence,
held-out results and provenance externally with hashes and one current entry.
Read-back checks preserve IDs, graph access and clinical joins. See
[data boundaries](data-boundaries.md) and [runtime](../runtime.md).

Sources: [scVI](https://docs.scvi-tools.org/en/stable/user_guide/models/scvi.html),
[scPoli implementation](https://github.com/theislab/scarches/blob/master/scarches/models/scpoli/scpoli_model.py),
[scPoli paper](https://www.nature.com/articles/s41592-023-02035-2).
