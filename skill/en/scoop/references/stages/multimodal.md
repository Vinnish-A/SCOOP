# Multimodal inventory, preparation and alignment

First organize exact IDs and provenance-linked dataset/sample/patient/library,
assay/platform and condition mappings under the
[metadata and batch contract](../docs/design-principles.md#metadata-and-batch-contract).
Assess batch handling for each measurement and analysis separately; RNA, ATAC,
ADT/protein and imaging do not share a correction by default. Unknown fields limit
only the methods or conclusions that require them.

Use for project-wide data preparation and for supporting modalities alongside
single-cell/spatial RNA. These are analytical contracts and validated adapter
routes; they do not imply that SCOOP has a native CLI module for every assay.
Read only the routes present in the project. Preserve the original measurements
and use a representative actual backend probe before expensive new processing.

## Inventory and reuse first

Search the live manifest **and** original raw/prepared resources, archive indexes,
source-study supplements and prior preprocessing receipts. A task-specific export
may omit whole assays, samples or cell populations. Recover omitted source assets
without automatically making them eligible for the current biological question.
Inspect archive members and object assays: a file called protein/counts/multiome
may contain TIFF images, normalized expression or separate unpaired libraries.
Check actual file/container headers before declaring readability: an `.xlsx`
download can be an HTML error page while a valid duplicate exists elsewhere.
Retain failed-source evidence and use the verified original numeric object.

For each source record the native object and assay/layer; species, tissue and
condition; raw/normalized/intensity/processed-count semantics; full feature and
observation axes; donor, specimen, capture, section, ROI/FOV and coordinate frame;
pairing and overlap evidence; existing QC/models; current eligibility and gaps.
Keep clinical genotype context distinct from actual cell-linked variant calls.
Gene-by-sample mutation-category tables may lack variant positions, alleles and
coverage; nonempty category records are not variant counts and blank is not
automatically wild type. Build separate exact sample-column maps for each assay,
even when their sample sets match. Preserve native row IDs for duplicated or
spreadsheet-date-like gene names; flag possible annotation damage and recover
authoritative original identifiers before merging or renaming them.
Source-view counts, repeated specimens and multiple modalities are not additive
independent patients or cells.
Identify each view by container **and native object/assay/layer**, not by physical
path alone. A single RData/workbook can hold RNA counts, VST, protein and other
matrices; retain distinct logical keys and axes instead of overwriting one view
when another shares its file.
Conversely, compatible old/new assay names can refer to the exact same native
object. Reconcile such aliases only after exact source/object/axis validation,
preserve the old-to-current resource-key map, and avoid registering them twice.
Directory indices, invalid downloads and their valid measurement children must
have explicit distinct roles; file existence alone does not establish an assay.

Inspect feature types inside multi-assay containers before QC or normalization.
A 10x H5 can store Gene Expression, Peaks and/or Antibody Capture on one barcode
axis even when the acquisition ledger calls it single-cell counts. Register each
assay separately, preserve the shared complete barcode/capture bridge, and use
only the relevant assay's features for depth, HVG and likelihoods. Never count
peaks or antibodies as RNA genes, or assume a separately exported RNA view adds
new cells. Determine whether recovered RNA is already represented in another
source view before expanding the cohort.

Use an executable per-resource ledger, for example:

`resource_id, modality, native_path, measurement, observation_unit,
subject_sample_capture_keys, feature_or_genome_contract, pairing_evidence,
valid_existing_artifact, reuse_status, missing_stage, next_action, receipt`.

Statuses distinguish verified reuse, executable missing processing, unresolved
source evidence and downstream dependencies. Verify input axes/signatures,
parameters, implementation and output hashes before inheriting a checkpoint.
Do not repeat valid complete-library QC or discard valid numeric fits because an
entry link moved. Changed scope invalidates only affected descendants. Per-cell
normalization caches can be subset after QC; a representation fitted on changed
cells/features needs a separate reuse decision. See [workflow receipts](../workflow-contract.md).

## Chromatin: ATAC, multiome, ChIP and CUT&Tag

- Record assembly, chromosome names, coordinate convention, interval IDs, feature
  ordering and peak-matrix semantics. BED zero-based half-open and GRanges
  one-based closed coordinates require an explicit conversion; width arithmetic
  alone does not establish assembly or provenance. Retain peak-ID crosswalks.
  Inspect the actual record layout; a file named BED may contain only literal
  `chr:start-end` IDs and does not establish an endpoint convention by extension.
- For same-nucleus RNA–ATAC, join complete native library+barcode IDs and verify
  one-to-one axes and donor/site/library attributes. Keep unmatched RNA-only and
  ATAC-only observations. Separate GEX/ATAC libraries from the same tissue do not
  establish same-cell pairing; statistical alignment remains an estimate.
- Reuse native fragment/library QC. When fragments and correct references are
  available, assess fragment depth, TSS enrichment, FRiP, nucleosome signal,
  blacklist and doublet evidence with an assay-appropriate implementation.
  A peak-only matrix cannot recreate absent fragment QC. Do not label missing
  QC as passed. Review processed/binarized/capped matrices before count models.
- For an ATAC representation, TF–IDF/LSI is a practical fast route; inspect depth
  association of components instead of automatically trusting or deleting one.
  For condition comparisons aggregate eligible native measurements by biological
  replicate × reviewed state, keeping source-specific peak universes and count
  semantics. Do not substitute per-cell P values for donor-level evidence.
- Motif analysis requires the correct sequence build, motif collection and
  GC/accessibility-aware background. Similar ETS motifs do not uniquely identify
  PU.1. Accessibility, motif deviation, peak–gene association and measured binding
  answer different questions. A matched peak ID does not establish regulation.
- For ChIP/CUT&Tag, retain antibody/target, assay, biological replicates, controls,
  peaks/signal/raw-read levels and assembly. Reuse supported peak/IDR calls or
  process reads with a documented assay pipeline when available. Record unmapped
  and multi-mapped liftovers. Do not turn promoter proximity into a confirmed
  target or binding in another species/cell model into binding in current cells.
  Use the source build's own chain, preserve native and target coordinate
  contracts, and inspect strand/width changes. A lifted bounding span crossing
  target insertions or gaps does not imply binding throughout that span; retain
  stretch flags and inspect aligned blocks before overlap interpretation.

Primary implementation references: [Signac ATAC workflow](https://stuartlab.org/signac/articles/pbmc_vignette),
[10x ATAC output/QC](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest/analysis/outputs/web-summary).

## Cell protein, CITE-seq, IMC and tissue proteomics

- Preserve native epitopes, antibody clones/isotypes, protein groups and explicit
  epitope-to-gene mappings. A family/isoform-reactive antibody is not a specific
  RNA gene. An absent target is unmeasured, not zero expression.
- CITE RNA–ADT pairing uses complete native capture-qualified IDs. Distinguish
  identity correspondence from verified gene-value equivalence to a separate RNA
  export. Keep genuinely unmatched ADT cells; do not rescue through a similar
  barcode in another library. Already connected CITE cells are not extra RNA cells.
- Determine whether ADT values are raw counts, CLR/dsb or other provider values.
  Normalize raw ADT using a suitable method such as CLR or, with appropriate raw
  background/control inputs, dsb. Record variant and inputs; do not apply raw-ADT
  normalization again to provider-normalized values. RNA UMI QC/likelihoods do
  not apply to compensated IMC intensity or mass-spectrometry abundance.
- IMC/multiplex images require channel/marker identity, compensation status,
  segmentation, image and coordinate frame, scale, batch and ROI/sample/patient
  nesting. Check native values, controls, spillover and segmentation evidence.
  Do not compensate twice. Explicitly transformed intensity may support phenotype
  review, not RNA-count interpretation or fabricated same-cell RNA pairing.
- Inspect every native TIFF page/series and preserve page indices. A GeoMx
  protein assay may export only morphology images while quantitative targets and
  controls were counted separately by barcodes. Page count is not protein-panel
  coverage; an unlabelled page cannot be assigned by expected marker order.
  Require measured ROI counts/controls and ROI area for the applicable numeric
  route, and explicit pixel scale/section registration for spatial image claims.
- Bulk proteomics requires protein/peptide grouping, batch/normalization and
  missingness assessment. Preserve NA and detection/quantification limits;
  distinguish plausible random from intensity-dependent missingness before
  choosing imputation or sensitivity analyses. Do not silently zero-fill, convert
  log intensity to RNA counts, or infer protein activity from RNA concordance.

References: [Seurat paired RNA/ADT workflow](https://satijalab.org/seurat/articles/multimodal_vignette),
[dsb paper](https://www.nature.com/articles/s41467-022-29356-8),
[DEP missingness workflow](https://www.bioconductor.org/packages/devel/bioc/vignettes/DEP/inst/doc/MissingValues.html).

## Bulk/regional RNA, perturbations and clinical outcomes

Determine count/TPM/FPKM/CPM/log/intensity semantics from provenance and native
values, not filenames. Preserve native feature identifiers and duplicate mapping
policy. Exact sample IDs and expression-column order must match clinical tables;
same specimen, same patient, same cell and independent cohort are different links.
Keep complete source tables and explicit eligibility/exclusion columns.

Raw sequencing counts permit count-based normalization and replicate-aware DE.
Normalized expression needs an appropriate continuous-data method and its limits;
never round TPM/CPM into pseudo-counts. FastDE's documented pseudobulk route remains
a fast option for compatible counts; use R edgeR/DESeq2 or limma as appropriate
reference implementations, not a universal automatic rerun. Clinical composition
estimates await a reviewed compatible reference; bulk deconvolution is distinct
from [spatial deconvolution](spatial.md).

Fractional values alone do not prove normalization: Salmon/RSEM/tximport gene
estimated counts can be fractional. Preserve their quantifier, gene IDs/version,
aggregation and length-offset provenance, and use a compatible count-analysis
route. Distinguish them from normalized/RLD values and integer single-cell UMIs;
do not reject or round every fractional bulk matrix by a universal rule.

For perturbations retain the real species, model, target, reagent/dose/time,
control, batch, library and experimental/biological replicate structure. Check
intervention efficacy, count depth, design rank and contrasts. Do not invent donor
pairing or treat technical libraries as donors. Preserve native IDs; cross-species
programme comparison requires a recorded orthology mapping with ambiguous/failed
members, not uppercase symbol conversion. Correct over the tested gene/programme
family. A cell-line perturbation supports its measured response; it does not
establish the full tissue interaction, secretion or direct binding mechanism.

For clinical outcomes retain endpoint, time origin/unit, censoring/event codes,
follow-up and treatment chronology. A recurrence outcome is not the specimen's
collection stage. Resolve repeated sampling and paired pre/post designs before
patient-level inference; ROIs, cells and bins do not supply independent survival
replicates. Unknown patient linkage need not block descriptive source-local work.
See [DE](differential-expression.md) and [outcomes](outcome-association.md).

Reference: [DESeq2 count/design contract](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html),
[limma user guide](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf).

## Execution and acceptance

Process missing work with sufficient source evidence while other branches wait
for labels or reference review. Preserve fast validated methods and bound CPU,
RAM and I/O; distinguish discovery/import checks from actual numeric execution.
Use selective reads and caches rather than repeatedly loading huge matrices.

Deliver the resource/reuse ledger; source and processed entry points; exact
join/feature/coordinate mappings; actual QC and model receipts; measured-feature
coverage; and residual gaps. Acceptance separates native identity/axis validity,
measurement/QC validity, biological annotation and fitness for a specific claim.
Do not mark the entire project aligned because one layer passed, or require all
assays to share one feature matrix/embedding. New reviewed labels invalidate
label-derived aggregates, references and outcomes, not unchanged native bridges.

## Omitted sources and unfiltered barcode matrices

Reconcile the current resource ledger against the last complete inventory and
the original acquisition receipts, including reference/model-organism sources
outside the active cohort. Audit actual files within source directories; a folder
or file count is not an assay count. Record original companion files and known
source-view overlap without promoting registry tables, model caches, label-based
subsets or rendered outputs to independent assays. A previously overlooked assay
does not automatically join a frozen model input.

Unfiltered barcode columns and positive-depth barcodes are not called cells.
Preserve the complete available barcode/background and feature axes, author group
labels and the original count package. Reuse complete stored-value/CRC checks
when current source hashes and axis/parameter boundaries agree; do not rescan
valid huge matrices just because the resource was omitted from a later ledger.
Cell calling is separate from numeric matrix QC and subsequent doublet analysis.
Verify chemistry, capture boundaries and available original calling evidence
before selecting a cell-calling route or applying droplet Scrublet. Directory
names do not establish precise drug combinations, animal counts or biological
replicates. Keep unresolved design fields unresolved and model species separate.

Infer MatrixMarket orientation by matching the actual barcode/observation and
feature axes, not by a presumed genes-by-cells convention. A real-valued header
does not by itself disqualify integer counts: check the actual full values and
source measurement provenance. Preserve duplicate literal symbols and missing
original gene accessions. Retain author-reported QC metrics beside newly computed
available-axis metrics and quantify discrepancies; exact cell IDs do not prove
unchanged gene/count universes. A fraction-unspliced metadata column or a SCENIC
loom is not a spliced/unspliced velocity measurement without the actual layers.

## Microarrays, methylation and external cohort validation

Keep native probes and full platform annotations, including controls and
one-to-many gene mappings. A gene-collapsed export cannot recreate original probe
measurements. Distinguish mRNA, microRNA and targeted panels, and preserve the
provider preprocessing scale. Do not repeat normalization merely because another
file format exists. Record provider imputation/batch removal for protein views;
these are different measurements from raw protein groups or RNA counts.

CpG methylation follows its own route: retain probe IDs, platform/build, sample
keys, controls and missingness. Confirm beta/M/other scale before transformation;
a [0,1] range alone does not prove beta provenance. Missing CpGs and measured zero
are different. A gene-level summary does not replace CpG measurements or RNA
expression. Processed matrices cannot supply absent raw IDAT/control QC. See
[minfi's native input and QC contracts](https://bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.html)
and [limma's omics scope](https://bioconductor.org/packages/release/bioc/html/limma.html).

Build exact clinical bridges at each assay's native sample/aliquot level. Missing
tokens such as literal `NA` stay in the original export but cannot be identity
keys. A participant-level join does not establish the same aliquot, time point
or exposure; preserve the [TCGA barcode hierarchy](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/).
Generated P1/P2 labels on different platforms cannot establish pairing.

Validate survival time units, event/censor definitions, endpoint and time origin
from source evidence. An old conversion script or forced all-event vector does
not establish them. Keep clinical IDH calls separate from transcriptional IDH
signatures. Compare original/derived assay views and cohort overlap before
independent validation; assay columns are not automatically independent patients.
External bulk evidence does not supervise integrated single-cell identities.
