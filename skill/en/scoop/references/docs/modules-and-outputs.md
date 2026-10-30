# 02. Modules and Deliverables

## qc

Inputs: `adata.layers['counts']`, `obs['sample_id']`, and gene symbols.

Default tools: Scanpy QC metrics, per-sample robust MAD/quantile thresholds, and `scanpy.pp.scrublet`.

Outputs:

- H5AD `obs`: QC metrics, QC flags, doublet scores, and doublet calls;
- External tables: `qc_thresholds_by_sample.tsv`, `scrublet_scores.parquet`;
- Figures: QC violins/scatters and Scrublet histograms;
- Logs: thresholds and removal fractions for each sample.

Do not: remove cells using a global `mt > 10%` rule; remove ribosomal-high or proliferation-high cells by default; treat Scrublet scores as the sole doublet evidence.

## core

Input: post-QC H5AD.

Default tool: FastCore. The native planner selects `fastcore_cpu`, `fastcore_mixed`, or `fastcore_oom`, with `scanpy_legacy` as the built-in fallback (Scanpy normalization/HVG/PCA/kNN/UMAP/Leiden plus Harmony 2.0). This is not an exclusive restriction on project adapters; if an entry point is unusable, bypass it according to the [first principles](design-principles.md) and record the actual backend, without claiming a native run.

Outputs:

- H5AD `layers['log1p_norm']`;
- H5AD `obsm`: biology PCA, identity PCA, Harmony PCA, and UMAP;
- H5AD `obs`: final clusters;
- External tables: HVG ranks, PC covariate associations, Leiden sweeps, and cluster stability.

Leiden uses a coarse-to-fine search within `0.25-1.5` by default: first identify stable candidates with a small set of resolutions/seeds, then perform a fine sweep around candidates using the full seed set. This avoids the cost of a full dense grid without treating one fixed resolution as the default answer.

Do not: retain clusters from every resolution/seed in `obs`; use UMAP distances for annotation; perform batch correction on condition.

Additional FastCore outputs:

- External JSON: `core/fastcore/fastcore_manifest.json`, `core/fastcore/core_quality.json`;
- Run logs: selected backend, fallback backend, whether fallback was used, and the quality-file path;
- Registry: register the FastCore manifest and quality report in `adata.uns['file_registry']['artifacts']`.

## programs

Input: log-normalized expression or a non-negative matrix derived from counts.

Default tool: `sklearn_multistart_nmf`, using independent random initializations, parallel K/seed fitting, component consensus, and usage refitting with fixed components. The standalone `fastcnmf` CLI is a separate implementation providing cNMF-compatible full consensus; do not conflate them. OmicVerse `single.cNMF` is an explicitly optional validation route.

Outputs:

- H5AD `obsm['X_nmf_usage']`, if programme and cell counts are manageable;
- External tables: NMF gene weights, K sweep, and programme summary;
- H5AD `obs`: optional dominant programme and programme entropy.

Do not: run every method in an NMF zoo by default; treat a single NMF run as a stable gene programme.

## annotation

Inputs: final clusters, marker table, NMF usage, kNN graph, and reference mapping.

Default tools: `04a_annotation_evidence.py` exports marker, programme, CNV, reference, and curated skill evidence; Scanpy markers are the default marker route, FastDE sparse COSG may provide fast auxiliary markers, and OmicVerse `single.find_markers(method='cosg')` is only optional compatibility validation. `04c_annotation_export_for_agent.py` only packages evidence into a structured template editable by a subagent/human; it does not generate final labels. `04d_annotation_commit.py` only validates and commits accepted structured decisions.

Outputs:

- H5AD `obs`: `cell_type_lvl1/2/3`, `cell_state`, and confidence;
- External tables: cluster markers, annotation evidence bundle, decision template, validation table, and merge/split log.

Do not: let marker-statistics scripts, reference models, or Python rules automatically decide cell types; let a subagent write directly to H5AD; rely on a single marker database, reference model, GPT, or CellVote alone to determine labels. Subprocesses may organize evidence and propose candidates; the main process must personally read the actual evidence and adjudicate final cell types, then pass them through the schema, validator, and committer.

### 04b_tumor_fastcnvpy

Trigger: the sample is tumor tissue and first-level `Major Lineage` annotation has been completed.

Default logic:

At runtime, require `annotation_confidence=high` or an explicit Boolean `reference_eligible_key`, combined with an identity allowlist and final_use to select references. By default, retain sparse matrices, compute in cell chunks, and disable N² cell clustering; explicitly enable `--compute-clusters` if needed. `chr_arm` accepts p/q or chromosome-matched forms such as 19.q, standardized to p/q.


- First distinguish parenchymal and non-parenchymal cells using `cell_type_lvl1`;
- Put a reviewed allowlist of high-confidence non-parenchymal identities into `fastcnv_reference_pool = normal_nonparenchymal`, balancing donors/identities and excluding mixed/unresolved/low-quality cells; do not directly include all non-epithelial cells as normal references;
- If normal parenchymal cells can be separated from tumor/candidate cells, place them in `normal_parenchymal`;
- If normal parenchymal and tumor cells are mixed, still run FastCNVpy, using a pooled non-parenchymal reference to provide CNV evidence for candidate parenchymal cells;
- FastCNVpy constructs a pooled reference on the merged H5AD, but splits by `sample_id` and computes CNV per sample to avoid misclassification caused by reference shifts in individual samples.

Outputs:

- H5AD `obs`: `fastcnv_reference_pool`, `fastcnv_cnv_fraction`, `fastcnv_normal_threshold`, `fastcnv_tumor_evidence`;
- External tables: pooled manifest, pooled cell metadata, pooled chromosome-arm CNV, genomic windows, and per-sample FastCNVpy results;
- Logs: reference composition, whether normal parenchymal cells were included, candidate cell counts, parallelism parameters, and dense/sparse H5AD mode.

Select parallelism and dense/sparse paths from an actual pilot; do not make old benchmarks the defaults for large objects. Control total worker memory/threads and avoid whole-object copies and N² cell clustering; disable nonessential clustering when necessary, process chunks against a fixed reviewed reference, and save checkpoints. Keep continuous/arm-level CNV evidence separate from final malignancy judgments; low-CNV cells may still be malignant.

## spatial

Inputs: spatial H5AD view, reference H5AD, coordinates, and spatial metadata.

Default tools: RCTD-py for deconvolution, Squidpy/spatial graph, and optionally BANKSY for domains.

Outputs:

- H5AD `obs`: spatial unit metadata, deconvolution summary, and domain labels;
- H5AD `obsm['spatial']`;
- External tables: RCTD weights, mode sensitivity, and spatial domain table;
- SpatialData/Zarr or GeoParquet: images, polygons, and transcript points.

Default rule: use RCTD-py `full` for low-resolution spots/ROIs; when each unit is expected to contain only a few dominant cell types, also run `multi` as a sensitivity analysis.

## ccc

Inputs: annotated H5AD, cell-type groups, and LR resource.

Default tool: FastCCC. Validate complex ligands/receptors with CellPhoneDB v5 or LIANA, which may be called through OmicVerse wrappers.

Outputs:

- External tables: FastCCC all/significant/differential results and complex-sensitive validation;
- H5AD `uns`: method summary, paths, hashes, and number of significant results;
- Optional H5AD `obs`: cell-type-level sender/receiver summaries.

Do not: write all FastCCC results into `adata.uns`; treat CCC as proof of mechanism; ignore subunit expression in multimeric ligands/receptors.

## differential-expression

Inputs: raw counts, cell type/state labels, and donor/sample/condition metadata.

Default tool: FastDE pseudobulk DESeq2-like negative-binomial Wald testing. Retain R DESeq2 and edgeR quasi-likelihood as reference validation. Sample-level cell abundance/state-proportion comparisons are handled separately by `fastde abundance`.

Outputs:

- External tables: pseudobulk counts, metadata, FastDE design matrix, FastDE DE table, size factors, dispersions, benchmark/R reference validation table, and GSEA;
- External tables: `fastde abundance` sample-by-celltype counts/proportions, mode-specific results, predictions, metrics, and manifest;
- H5AD: only result indexes and necessary summaries.

Do not: use Harmony/scVI/integrated expression for DE; use per-cell p values to establish condition-level differences; interpret marker-gene tests as condition DE; interpret abundance as gene-level DE; interpret the current FastDE as a complete replication of `lfcShrink`.

FastDE condition DE output columns retain the DESeq2 style: `gene`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`, `dispersion`, with additional `dispGeneEst`, `dispFit`, `dispMAP`, and `dispOutlier` columns. In the current full-dispersion benchmark with 30k genes / 12 pseudobulk samples, FastDE takes approximately `18.55s` and the R DESeq2 reference approximately `19.68s`; log2FC Spearman correlation is approximately `0.99999`, and top100 overlap is `1.00`.

FastDE abundance output columns vary by mode: survival writes `hazard_ratio_approx` and concordance; binary writes `odds_ratio_approx`, positive/negative mean proportions, and classification metrics; multiclass writes class/reference contrasts. All abundance inference is sample-level.
