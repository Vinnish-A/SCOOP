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
    "fastcopykat.cnv_prediction": EngineSpec(
        engine_id="fastcopykat.cnv_prediction",
        task_type="cnv_prediction",
        version="0.1.0",
        input_schema="h5ad.raw_counts.v1",
        output_schema="fastcopykat.cnv_prediction.v1",
        consumes=("h5ad", "counts_layer", "sample_key"),
        produces=("cnv_prediction_table", "cnv_manifest"),
        default_cli=("fastcopykat", "predict"),
        writes_h5ad_fields=("fastcopykat_prediction",),
        writes_external_artifacts=("fastcopykat_predictions.tsv", "fastcopykat_manifest.json"),
        quality_gates=("min_cells", "reference_quality"),
    ),
}


def list_engines() -> list[EngineSpec]:
    return [_ENGINES[key] for key in sorted(_ENGINES)]


def get_engine(engine_id: str) -> EngineSpec:
    try:
        return _ENGINES[engine_id]
    except KeyError as exc:
        raise KeyError(f"unknown SCOOP Fast engine: {engine_id}") from exc
