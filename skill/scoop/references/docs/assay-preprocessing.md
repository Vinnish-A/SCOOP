# Assay-specific preprocessing

Decision updated 2026-10-16: retain the ordinary scRNA workflow. Imaging spatial
assays need a simple starting range and a QC checklist, not a universal hard
threshold or an extra approval ceremony. Select `run.assay` explicitly.

| Assay | Practical starting range | SCOOP permissive starting values |
| --- | --- | --- |
| Xenium | 10–50 target transcripts and 5–20 detected target genes per cell | 10 transcripts, 5 genes |
| CosMx | 20–50 target transcripts and 10–30 detected target genes per cell | 20 transcripts, 10 genes |
| Visium / legacy spot ST | Review counts/genes and tissue membership as spots, not single cells | 100 counts, 50 genes |

These ranges are engineering starting points, not validated biological cutoffs.
Use `qc.spatial.min_counts` and `min_genes` to override them. Small custom panels,
cell size, section thickness, and the local cell population change their meaning.
Do not apply the same gene thresholds to a 50-gene and a 5K panel. Review the
retained fraction and its tissue/FOV distribution before tightening thresholds,
especially when studying rare immune cells such as DCs.

Check these metrics:

- Target transcript counts and detected genes; keep controls out of both totals.
- Available negative probe/codeword or other non-gene-feature counts. Missing
  controls mean unassessed background, not zero background.
- Cell/nucleus area, transcript density, geometry and segmentation boundaries;
  investigate mixtures with morphology/marker evidence, not default Scrublet.
- Finite coordinates, tissue membership, FOV/sample shifts, and spatial patterns
  of excluded observations. Available FOV metadata (`qc.spatial.fov_key`, default
  `fov`) produces a separate sample/FOV summary. Optional `qc.spatial.cell_area_key`, `in_tissue_key`
  and `source_pass_key` apply the corresponding source metadata.
- Molecular decoding quality when supplied upstream. A processed expression
  matrix cannot reconstruct transcript-level quality scores or segmentation.

No mitochondrial/ribosomal filter or droplet-doublet detection is applied to the
imaging spatial route. Panel mitochondrial percentages are descriptive only;
missing targets are recorded as missing rather than zero. Native counts and
source observations are retained by QC; the spatial core uses `final_use`.

## Features and normalization

The spatial core is a simple baseline, separate from RNA FastCore/Scanpy. It
starts with expressed target genes and library normalization/log1p; it does not
calculate whole-transcriptome stress/cycle scores. `core.spatial.features` can
be `all_expressed` or `seurat` (dispersion HVGs on log expression). Larger 5K/6K
panels warrant comparing these routes; targeted assays do not imply "never HVG".
`core.spatial.normalization: log1p` omits library-size normalization for sensitivity
analysis. Tissue-domain detection and cell-state clustering have different aims;
no normalization is claimed universally optimal. These settings never replace
raw counts for deconvolution or count models. Batch correction defaults to none
in the spatial route; it requires an explicit spatial setting to enable Harmony.

## Small scRNA lineage / LOESS failure

Seurat v3 on integer counts remains the ordinary RNA primary method. Calculate
HVGs in a lightweight feature-metadata object so inherited mixed `var` columns
cannot break Scanpy's batch concatenation. Keep all original counts and genes;
only zero/unselected features are excluded from PCA selection.

Known LOESS singularity/extrapolation failures or missing `scikit-misc` permit an
explicitly recorded Seurat dispersion fallback on log-normalized expression.
Unrelated exceptions propagate. Optional `core.hvg_retry_spans` permits trying a
larger span; the default is no extra retries, since increasing span can change
which genes are selected. `core.hvg_fallback: error` disables the fallback.
Record actual method, failures, span, input representation and per-batch sizes
in `uns.hvg_selection` and the external HVG table. Do not drop donors, merge
biological groups, or declare low-quality cells merely to suppress LOESS errors.

## Evidence consulted

- [MacDonald et al., Immunity 2024](https://doi.org/10.1016/j.immuni.2024.11.004)
  and [author methods/code](https://github.com/DrLucyMac/spatial_synoDC): the
  corresponding RA CosMx methods use 30 counts and 20 genes after segmentation.
  The accessible preprint reports these thresholds explicitly; retain that
  version distinction when comparing sample/FOV counts across releases.
- [Salas et al., Nature Methods 2025](https://www.nature.com/articles/s41592-025-02617-2):
  evaluates preprocessing of small Xenium panels; normalization, scaling and
  feature-selection choices affect clustering. Its parameters are not universal
  defaults for 5K panels or RA.
- [Systematic platform benchmark, Nature Communications 2025](https://www.nature.com/articles/s41467-025-64990-y):
  uses permissive transcript minima of 10 for Xenium and 20 for CosMx.
- [IBD atlas, Nature Communications 2026](https://www.nature.com/articles/s41467-026-72482-w):
  uses platform/panel-specific QC and 3,000 HVGs for the 5K/6K panels, illustrating
  why a blanket ban on spatial HVGs is inappropriate.
- [Bioconductor OSTA, imaging QC](https://bioconductor.org/books/3.22/OSTA/pages/img-quality-control.html):
  controls, morphology, FOV/tissue context, and exclusion bias matter. Mitochondrial
  and ribosomal panel coverage can be inadequate for standard RNA QC.
- [10x Xenium output specification](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/1.9/analysis/outputs/xoa-output-understanding-outputs):
  the default cell-feature matrix contains Q20-qualified transcripts and explicit
  feature types; these controls must not enter biological expression analyses.
- [Bhuva et al., Genome Biology 2024](https://doi.org/10.1186/s13059-024-03241-7):
  demonstrates tissue-domain information in library sizes and motivates retaining
  unnormalized/log-only sensitivity analysis for domain detection.
- [Scanpy HVG contract](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html)
  and [maintainer discussion #1504](https://github.com/scverse/scanpy/issues/1504):
  v3 consumes counts; dispersion methods consume log data. Small batches and
  low-expression distributions can cause LOESS failures; filtering or a larger
  span may help, but neither establishes a universal repair. This is developer
  evidence, not a published benchmark establishing optimal recovery parameters.
