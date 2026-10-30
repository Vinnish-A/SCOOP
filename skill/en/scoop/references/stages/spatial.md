# Spatial views, graphs and RCTD-py

Apply the [metadata and batch contract](../docs/design-principles.md#metadata-and-batch-contract)
to the reference and each native spatial assay, platform, specimen/section and
FOV. Declare the method's actual calibration/batch interface and assessed scope;
RNA latent mixing does not calibrate spatial counts or establish registered views.
Keep native likelihood inputs and exposure, and preserve separate sections and
unresolved pairing where the method cannot support joint correction.

Inputs:

- Spatial H5AD view.
- Annotated single-cell reference H5AD.
- Spatial coordinates and optional SpatialData sidecar.

Script:

```bash
scoop spatial --config runs/<run_id>/config/run.yaml --dry-run
scoop spatial --config runs/<run_id>/config/run.yaml
```

## Reference granularity: all cell subgroups

Spatial deconvolution must resolve **all reviewed cell subgroups across every
lineage in the tissue**, rather than splitting only the target cells while
collapsing the remaining reference into T, B, myeloid or stromal classes. Use the
finest evidence-reviewed identities available for each compartment. Fit this
complete fine reference first; sums of its subgroup weights can support broad
views or composition adjustment afterward.

Before fitting, export a reference coverage table linking each subgroup to its
parent lineage, reviewed identity source, cell and donor counts, inclusion status
and any unresolved gap. Preserve rare, mixed or genuinely unresolved populations
under explicit labels. If a required subgroup is missing or cannot be separated,
record that coverage or identifiability gap and prepare the missing reference or
annotation evidence; do not claim a broad-only fit resolves all subgroups. A
coarse reference may be retained as a declared pilot or sensitivity analysis.
Do not create subtypes through automatic label transfer to satisfy the table.

Inspect subgroup profiles and marker separability, especially closely related
states. Overlapping profiles call for a recorded merged-state sensitivity
analysis alongside the fine fit, not silent removal of inconvenient subgroups.

## Execution

Default rule:

- Low-resolution spot/ROI/large bin with compatible native RNA counts: RCTD-py `full` mode.
- Sparse or near-single-cell units: run `multi` or `doublet` only when justified.
- Check the specific backend with a representative real pilot, not only a
  `spatial: ok` flag. Use a recorded project adapter if the CLI/backend is unusable.

Classify the measurement before choosing the likelihood. A GeoMx ROI table may
already contain Q3-normalized RNA or area-normalized protein; its image archive
is a different measurement view. Recover native counts/controls/exposure if
available. Otherwise preserve the provider scale and validate a method compatible
with that input (and its identifiability limits); do not round normalized values
or pass them unchanged to count-likelihood RCTD/cell2location. Protein ROI values
are not a spatial RNA matrix. Repeated ROI labels need native row/segment keys,
not automatic merging; specimen correspondence does not establish registered ROIs.

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

## Version, reference and speed validation

The runtime pin above describes the installed backend, not equivalence to every
upstream release. [Upstream's R comparison](https://rctd-py.readthedocs.io/en/latest/r_spacexr.html)
records a pre-0.4.0 bulk-calibration gene-selection difference and changes the
sigma sampling default `N_fit` from100 to1000. Do not cite current upstream parity
as validation of an older pin or a project using a smaller iteration/sample budget.
Record actual package/build, genes used at each step, full native UMI exposure,
sigma, sampling seed/budget, iterations/convergence and output type coverage.
Check which settings the actual mode implementation consumes: an adapter config
or receipt can contain an ignored `max_iter`. Inspect/test the executed API path
before attributing speed or convergence to such an override.

Keep the fast RCTD route. When a release difference or changed reference matters,
run representative real sections with the installed backend, a staged candidate
and/or R spacexr using comparable mode/parameters. Preserve the shared runtime
until an upgrade is separately built/validated; an isolated project adapter can
record candidate results. A successful run, speed or dominant-type agreement
alone does not establish accurate fine mixtures. Inspect abundance concordance,
residual/held-out genes, rare/related-type identifiability and tissue/pathology
context. Use a heavier method such as cell2location on informative representative
sections when it answers a specific unresolved concern, rather than automatically
doubling every whole-cohort fit. Report measured time and memory on the same input.

Preserve complete native reference counts and their total depth **before** any
feature intersection or marker selection. A small marker panel may be an explicit
pilot, not proof of sufficient fine-type separability. Check every intended type
against types actually retained by the backend: `cell_min`/UMI filters can silently
drop rare populations. Preserve posterior/uncertainty information where available;
do not remove unresolved/rare types and renormalize the rest without reporting a
separate sensitivity analysis. Held-out validation genes must be absent from
fitting/feature selection and marker-derived ground truth is not independent truth.

Reference estimators legitimately calculate per-gene count means/profiles within
reviewed types. This differs from assigning identities by a marker-score mean;
do not reject valid reference profile estimation merely because it uses means.
Keep the complete reference immutable across sections: take a fresh per-section
gene intersection rather than shrinking the shared reference in a loop, which
would make the fitted genes depend on section order.

When learning from published studies, verify the actual spatial methods and code
separately from bulk deconvolution or module scoring. A paper's BayesPrism or
CIBERSORTx bulk analysis does not establish its ST deconvolution method. Borrow
validated methodological choices, not undocumented subtype deletions or fitted
results. Native segmented-cell imaging ordinarily follows its measured panel/QC
and morphology, not spot mixture inference by default.


## Paired count views and unresolved reference components

An annotated normalized-expression UID can use a separately published real count view only after an explicit native row/feature bridge and numerical gene-value audit. Preserve annotation UID and count-source UID/row/feature identifiers separately. Read actual integer counts from the count view; do not fabricate them by exponentiating or rounding logCPM. State whether full depth means the available published filtered axis or an independently recovered original library. Recheck required markers and separability when this view changes the common reference gene axis; missing genes are not zero-filled.

Qualify the actual residual UID subset after refinements/arbitration. Whole-parent purity and DEG do not prove the remaining subset has that identity. Broad RNA can be reliable while adequate reviewed fine evidence remains conflicting; manually reviewed native-count components may retain that unresolved identity without rewriting atlas labels. Low-information or assay-limited populations are evidence-pending coverage/identifiability gaps: actively seek available evidence and record unavailable measurements, rather than calling missing evidence a resolved biological conflict. Near-convex overlap with known profiles is an identifiability warning, not a label classifier. Report conditional RNA allocations and aggregate/merged sensitivities; numerical convergence does not make every component an independently measurable cell proportion.


## Validation scope and format adapters

Restore a raw MatrixMarket/Space Ranger source through an explicit native adapter
before sending it to an H5AD reader. Verify feature-by-barcode orientation,
barcodes, source feature IDs and literal symbols; independently reread integer
counts and full row depth. Preserve original rows and source tissue/QC flags.
Matching spot IDs alone does not complete count-derived QC or format support.

Use held-out native cells or sample units to stress related and rare reference
components. If a cached common-gene matrix omits actually measured genes, preserve
their original count mass separately and restore the known full exposure before
calibration; never add a fictional catch-all gene. Constructed mixtures measure
recovery of recorded source RNA origins, not independent biological identity or
cell-number truth. Separate global sample-unit holdouts from within-unit technical
cell holdouts, report types with no holdout, and avoid patient-level intervals
from repeated cell draws.

Compare errors at fine and reviewed parent levels, absent-component allocation,
depth and profile overlap. Stable fast/standard or backend agreement can coexist
with poor fine recovery; more epochs need not fix an inseparable reference. Keep
the fine fit and explicit unresolved components, export declared parent/merged
sensitivities, and withhold unsupported fine-abundance interpretation. Complete
cohort fitting and biological interpretation are distinct delivery states.


Carry explicit reviewed source exclusions into spatial input **before platform
calibration and fitting**. Merely retaining flags or filtering plotted weights
leaves excluded observations in the learned platform effects. Keep original rows
in the delivered table with not-fitted/NA status; preserve an unmasked run only as
a declared sensitivity. A missing unified `analysis_use` field does not prove QC
was absent: join actual count-QC sidecars by exact source row/ID and verify depth.
Keep author retention, count-derived review flags, coordinate validity, model
filters, convergence and biological approval as separate fields rather than
synthesizing an undocumented universal QC-pass flag.


For model comparisons, align the outcome scale as well as genes and locations.
A cell2location cell-number posterior is not an RCTD RNA fraction. Compare
model-expected component RNA on a declared shared gene set when needed; retain
background handling, reference estimator and actual training-cell differences as
joint-workflow sensitivities. A ratio of posterior component means is not a
credible interval for a normalized proportion.

Persist actual objective traces independently of console progress and verify they
are nonempty after resuming. Disabling a metrics logger can remove the evidence
needed for convergence review; an empty model history may also miss metrics from
a resumed trainer. Use actual trainer/CSV records if necessary and distinguish
the recorded continuation from earlier unrecorded epochs. A flat global objective
can coexist with changing fine allocations and does not establish biological
identifiability.

## Managed cached-reference and recovery route

The locked API adapter uses rctd-py 0.4.0. `spatial.rctd.preset` is `fast`
(N_fit=50, N_epoch=4) or `standard` (1000,8); explicit N_fit/N_epoch overrides
are recorded. The upstream full solver uses 50 iterations; `max_iter` is not
an effective override. Set actual thread/device budgets for concurrent runs.
`gene_selection: all_shared` is an explicit sensitivity using the same shared
gene axis for calibration and fitting; the default remains data-derived DE genes.
The adapter verifies the actual retained axis. A shared-signature control may
pass the same empirical native-rate profiles (with a declared common scale) to
an independent model without claiming an NB reference regression was performed.
Compare modeled component RNA on the same genes and locations, including each
method's platform/background treatment, rather than equating cell-number
posteriors with full-native RNA fractions.

For a large already-qualified reference, provide `reference_profiles` (NPZ with
`profiles`, `gene_names`, `cell_type_names`) and `reference_manifest` under
`spatial.rctd`. The JSON manifest uses `schema_version: scoop.reference-profiles.v1`,
`unit: mean_counts_per_full_native_exposure`, `labels_reviewed: true`,
`native_full_exposure_preserved: true`, `profiles_sha256`, per-type
`training_cells`, and explicit `coverage_complete`/`gaps`. Bind exact UID/QC,
source axes and review receipts in the project construction evidence. This
manifest records analyst qualification; naming a component is not proof of it.

The H5AD route requires reviewed annotation status or an explicit reviewed
`reference_mask_key`. Declared types dropped by backend thresholds cause an
error; rare types require qualification rather than implicit deletion.
Set `spatial_counts_layer: X` or `reference_counts_layer: X` only for explicit
native-count X inputs. A projected count view can declare `spatial_exposure_key`
for its verified full RNA depth. Projected reference H5ADs require cached
native-exposure profiles instead of a constructor that renormalizes their depth.

An optional `spatial_registry` Parquet joins every original row by exact
`source_row`/`source_axis_id`. Explicit source masks apply before calibration;
unfitted original rows remain in weights/H5AD with NA. Missing mask values remain
unassessed. Count-derived QC in a separate sidecar must be joined explicitly.

The adapter saves a hash-bound computation checkpoint before tables/H5AD export.
Retrying the same fit reuses it after input, code, parameters, package and lookup
checks; changing export serialization alone does not refit. Source/reference
changes require a new computation. Successful output still has limited validation
and unassessed fine-proportion interpretation until the relevant science is reviewed.

When comparing deconvolution methods, separate estimator, actual training UID,
feature selection, calibration and reference coverage changes. Use controlled
same-cell/shared-gene comparisons before attributing discordance to an algorithm.
