from __future__ import annotations

from .engine_spec import EngineSpec


_ENGINES: dict[str, EngineSpec] = {
    "fastde.markers": EngineSpec(
        engine_id="fastde.markers",
        task_type="marker_genes",
        version="0.1.0",
        input_schema="h5ad.expression_groups.v1",
        output_schema="fastde.markers.v1",
        consumes=("h5ad", "obs[groupby]", "expression_layer"),
        produces=("marker_table", "marker_manifest"),
        default_cli=("fastde", "markers"),
        writes_h5ad_fields=(),
        writes_external_artifacts=("markers_cosg.tsv", "markers_cosg_manifest.json"),
        quality_gates=("min_cells_per_group", "nonempty_markers"),
    ),
    "fastde.pseudobulk_deseq2": EngineSpec(
        engine_id="fastde.pseudobulk_deseq2",
        task_type="condition_de",
        version="0.1.0",
        input_schema="pseudobulk.counts_metadata.v1",
        output_schema="fastde.deseq2_wald.v1",
        consumes=("counts.tsv", "metadata.tsv", "condition_col"),
        produces=("de_table", "design_matrix", "size_factors", "dispersions"),
        default_cli=("fastde", "deseq2"),
        writes_h5ad_fields=(),
        writes_external_artifacts=("de_fastde_deseq2.tsv", "fastde_deseq2_manifest.json"),
        quality_gates=("min_samples_per_group", "dispersion_fit"),
    ),
    "fastde.abundance": EngineSpec(
        engine_id="fastde.abundance",
        task_type="differential_abundance",
        version="0.1.0",
        input_schema="sample_by_celltype_counts.v1",
        output_schema="fastde.abundance.v1",
        consumes=("h5ad.obs[sample_key]", "h5ad.obs[celltype_key]", "sample_metadata"),
        produces=("abundance_results", "abundance_predictions", "abundance_metrics", "abundance_manifest"),
        default_cli=("fastde", "abundance"),
        writes_h5ad_fields=(),
        writes_external_artifacts=("abundance_*_results.tsv", "abundance_*_predictions.tsv", "abundance_*_metrics.json", "abundance_manifest.json"),
        quality_gates=("min_samples", "min_events_for_survival", "convergence", "nondegenerate_labels"),
    ),
    "fastcnmf.programs": EngineSpec(
        engine_id="fastcnmf.programs",
        task_type="gene_programs",
        version="0.1.0",
        input_schema="h5ad.nonnegative_expression.v1",
        output_schema="fastcnmf.programs.v1",
        consumes=("h5ad", "expression_layer", "k_grid"),
        produces=("program_weights", "usage", "k_sweep", "manifest"),
        default_cli=("fastcnmf", "run"),
        writes_h5ad_fields=("obsm.X_nmf_usage",),
        writes_external_artifacts=("program_weights.tsv", "k_sweep.tsv", "fastcnmf_manifest.json"),
        quality_gates=("stability_threshold", "replicate_count"),
    ),
    "fastcnvpy.tumor_pooled": EngineSpec(
        engine_id="fastcnvpy.tumor_pooled",
        task_type="tumor_cnv_evidence",
        version="0.1.0",
        input_schema="h5ad.tumor_lineage_reference.v1",
        output_schema="fastcnvpy.pooled_result.v1",
        consumes=("h5ad", "gene_metadata", "sample_key", "major_lineage_key"),
        produces=("cnv_evidence_obs", "pooled_manifest", "sample_results"),
        default_cli=("python", "scripts/04b_tumor_fastcnvpy.py"),
        writes_h5ad_fields=("fastcnv_reference_pool", "fastcnv_cnv_fraction", "fastcnv_normal_threshold", "fastcnv_tumor_evidence"),
        writes_external_artifacts=("tumor_fastcnv_fastcnvpy_pooled_manifest.json",),
        quality_gates=("reference_pool_size", "candidate_cell_count"),
    ),
    "fastcore.preprocess": EngineSpec(
        engine_id="fastcore.preprocess",
        task_type="core_preprocess",
        version="0.1.0",
        input_schema="h5ad.raw_counts.v1",
        output_schema="fastcore.preprocess.v1",
        consumes=("h5ad", "counts_layer", "sample_key"),
        produces=("log1p_layer", "hvg_table", "program_scores"),
        default_cli=("fastcore", "preprocess"),
        writes_h5ad_fields=("layers.log1p_norm", "var.highly_variable_biology"),
        writes_external_artifacts=("hvg_biology.parquet", "fastcore_manifest.json"),
        quality_gates=("hvg_nonempty", "no_nan_scores"),
    ),
    "fastcore.core_pipeline": EngineSpec(
        engine_id="fastcore.core_pipeline",
        task_type="core_pipeline",
        version="0.1.0",
        input_schema="h5ad.qc.v1",
        output_schema="scoop.core_h5ad.v1",
        consumes=("h5ad", "core_config", "backend_plan"),
        produces=("adata_core.h5ad", "fastcore_manifest", "core_quality"),
        default_cli=("python", "scripts/02_core_analysis.py"),
        writes_h5ad_fields=("obsm.X_pca_biology", "obsm.X_pca_identity_prebatch", "obsm.X_umap_identity", "obs.cluster_identity"),
        writes_external_artifacts=("fastcore_manifest.json", "core_quality.json", "cluster_stability.parquet"),
        quality_gates=("pca_quality", "graph_quality", "cluster_stability"),
    ),
    "fastcore.graph_embed_cluster": EngineSpec(
        engine_id="fastcore.graph_embed_cluster",
        task_type="graph_embed_cluster",
        version="0.1.0",
        input_schema="h5ad.pca.v1",
        output_schema="fastcore.graph_embed_cluster.v1",
        consumes=("pca_representation", "neighbors_config", "leiden_config"),
        produces=("connectivities", "distances", "umap", "cluster_identity", "leiden_sweep"),
        default_cli=("fastcore", "graph"),
        writes_h5ad_fields=("obsp.connectivities_identity", "obsp.distances_identity", "obsm.X_umap_identity", "obs.cluster_identity"),
        writes_external_artifacts=("leiden_sweep.parquet", "cluster_stability.parquet"),
        quality_gates=("neighbor_overlap", "umap_valid", "leiden_ari"),
    ),
    "fastcore.quality_compare": EngineSpec(
        engine_id="fastcore.quality_compare",
        task_type="core_quality_compare",
        version="0.1.0",
        input_schema="fastcore.comparison_inputs.v1",
        output_schema="fastcore.quality_report.v1",
        consumes=("candidate_h5ad", "reference_h5ad", "quality_config"),
        produces=("core_quality", "pca_quality", "graph_quality", "cluster_stability"),
        default_cli=("fastcore", "quality"),
        writes_h5ad_fields=(),
        writes_external_artifacts=("core_quality.json", "pca_quality.tsv", "graph_quality.tsv", "cluster_stability.parquet"),
        quality_gates=("min_pc_subspace_cosine", "min_neighbor_overlap", "min_leiden_ari_vs_reference"),
    ),
}


def list_engines() -> list[EngineSpec]:
    return [_ENGINES[key] for key in sorted(_ENGINES)]


def get_engine(engine_id: str) -> EngineSpec:
    try:
        return _ENGINES[engine_id]
    except KeyError as exc:
        raise KeyError(f"unknown SCOOP Fast engine: {engine_id}") from exc
