from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np
import pandas as pd

from fastcnvpy import FastCNVConfig, PooledFastCNVResult, run_fastcnv_pooled_anndata
from fastcnvpy.io import write_pooled_outputs
from scsp_agent_sop.storage import write_table


DEFAULT_NONPARENCHYMAL_LINEAGES = {
    "B",
    "T",
    "NK",
    "Myeloid",
    "Mast",
    "Plasma",
    "Endothelial",
    "Fibroblast",
    "Stromal",
    "Pericyte",
    "Immune",
}


@dataclass(frozen=True)
class TumorFastCNVPlan:
    reference_key: str
    reference_labels: list[str]
    n_reference_cells: int
    n_normal_nonparenchymal: int
    n_normal_parenchymal: int
    n_candidate_cells: int
    should_run_fastcnv: bool
    reason: str
    nonparenchymal_lineages: list[str]
    dropped_reference_by_lineage: dict[str, int]
    measurement: dict = field(default_factory=lambda: {"status": "not_declared"})


def prepare_tumor_fastcnv_reference(
    adata,
    *,
    major_lineage_key: str = "cell_type_lvl1",
    sample_key: str = "sample_id",
    parenchymal_lineages: list[str] | tuple[str, ...] = ("Epithelial",),
    nonparenchymal_lineages: set[str] | None = None,
    normal_status_key: str | None = None,
    normal_status_labels: list[str] | tuple[str, ...] = ("normal", "non_tumor", "healthy"),
    reference_key: str = "fastcnv_reference_pool",
    reference_eligible_key: str | None = None,
    min_reference_cells: int = 20,
) -> TumorFastCNVPlan:
    """Build the pooled-reference labels used by the tumor FastCNV step."""

    if major_lineage_key not in adata.obs.columns:
        raise ValueError(f"major_lineage_key {major_lineage_key!r} is not present in obs")
    if sample_key not in adata.obs.columns:
        raise ValueError(f"sample_key {sample_key!r} is not present in obs")

    obs = adata.obs
    major = obs[major_lineage_key].astype(str)
    parenchymal = set(map(str, parenchymal_lineages))
    nonpar = set(DEFAULT_NONPARENCHYMAL_LINEAGES if nonparenchymal_lineages is None else nonparenchymal_lineages)
    if reference_eligible_key:
        if reference_eligible_key not in obs or obs[reference_eligible_key].dtype != bool:
            raise ValueError("reference_eligible_key must identify an explicit non-missing boolean review mask")
        eligible = obs[reference_eligible_key].copy()
    else:
        if "annotation_confidence" not in obs:
            raise ValueError("CNV reference requires reviewed annotation_confidence or reference_eligible_key")
        eligible = obs.annotation_confidence.astype(str).eq("high")
    if "annotation_status" in obs:
        # Structured annotation commits emit "accepted"; "labeled" is legacy.
        eligible &= obs.annotation_status.astype(str).isin({"accepted", "labeled"})
    if "final_use" in obs:
        eligible &= obs.final_use.eq(True)

    is_parenchymal = major.isin(parenchymal)
    is_normal_nonpar = major.isin(nonpar) & ~is_parenchymal & eligible
    is_normal_par = pd.Series(False, index=obs.index)
    if normal_status_key:
        if normal_status_key not in obs:
            raise ValueError(f"normal_status_key {normal_status_key!r} is not present in obs")
        is_normal_par = is_parenchymal & eligible & obs[normal_status_key].astype(str).isin({str(x) for x in normal_status_labels})

    ref = pd.Series("candidate_tumor_or_mixed", index=obs.index, dtype=object)
    ref.loc[is_normal_nonpar] = "normal_nonparenchymal"
    ref.loc[is_normal_par] = "normal_parenchymal"
    adata.obs[reference_key] = ref

    labels = ["normal_nonparenchymal"]
    if int(is_normal_par.sum()) > 0:
        labels.append("normal_parenchymal")
    n_ref = int(ref.isin(labels).sum())
    n_candidate = int((ref == "candidate_tumor_or_mixed").sum())
    should_run = bool(n_ref >= min_reference_cells and n_candidate > 0)
    if not should_run:
        reason = "insufficient pooled reference cells or no candidate parenchymal/tumor cells"
    elif int(is_normal_par.sum()) > 0:
        reason = "normal parenchymal and tumor/candidate cells are separable; FastCNV will provide orthogonal CNV evidence"
    else:
        reason = "normal parenchymal is mixed or unavailable; FastCNV will identify CNV-supported tumor cells using pooled nonparenchymal reference"
    # Reviewed nonparenchymal cells outside the vocabulary are excluded, not hidden.
    dropped = major[eligible & ~is_parenchymal & ~major.isin(nonpar)].value_counts()

    return TumorFastCNVPlan(
        reference_key=reference_key,
        reference_labels=labels,
        n_reference_cells=n_ref,
        n_normal_nonparenchymal=int((ref == "normal_nonparenchymal").sum()),
        n_normal_parenchymal=int((ref == "normal_parenchymal").sum()),
        n_candidate_cells=n_candidate,
        should_run_fastcnv=should_run,
        reason=reason,
        nonparenchymal_lineages=sorted(nonpar),
        dropped_reference_by_lineage={str(k): int(v) for k, v in dropped.items()},
    )


def run_tumor_fastcnv_workflow(
    adata,
    gene_metadata: pd.DataFrame,
    *,
    output_dir: Path,
    sample_key: str = "sample_id",
    major_lineage_key: str = "cell_type_lvl1",
    parenchymal_lineages: list[str] | tuple[str, ...] = ("Epithelial",),
    normal_status_key: str | None = None,
    normal_status_labels: list[str] | tuple[str, ...] = ("normal", "non_tumor", "healthy"),
    layer: str | None = None,
    h5ad_mode: str = "sparse",
    reference_eligible_key: str | None = None,
    nonparenchymal_lineages: set[str] | None = None,
    all_cells_q75_fallback: bool = False,
    compute_clusters: bool = False,
    cell_block_size: int = 10000,
    n_jobs: int = 1,
    config: FastCNVConfig | None = None,
    sample_name: str = "tumor_fastcnv",
    measurement=None,
) -> tuple[PooledFastCNVResult | None, TumorFastCNVPlan]:
    """Pooled-reference FastCNV evidence.

    The pooled engine shares one gene selection, scale factor and threshold
    across samples, so with ``measurement`` (a ``FeatureMeasurement``) only
    genes measured in every sample of the run enter it; a gene some sample
    did not quantify would otherwise be a structural zero read as a loss.
    The excluded genes are written next to the CNV outputs and counted in
    ``plan.measurement``.
    """
    if not adata.var_names.is_unique:
        raise ValueError(
            "CNV genomic positions require unique analysis feature names; "
            "resolve duplicate symbols through an explicit locus/ENSG mapping before fitting"
        )
    plan = prepare_tumor_fastcnv_reference(
        adata,
        major_lineage_key=major_lineage_key,
        sample_key=sample_key,
        parenchymal_lineages=parenchymal_lineages,
        normal_status_key=normal_status_key,
        normal_status_labels=normal_status_labels,
        nonparenchymal_lineages=nonparenchymal_lineages,
        reference_eligible_key=reference_eligible_key,
    )
    if measurement is not None:
        measurement.check(adata)
        samples = adata.obs[measurement.sample_key]
        measured = measurement.measured_in_all(samples)
        present = np.unique(measurement.sample_codes(samples))
        positioned = adata.var_names.isin(gene_metadata["hgnc_symbol"].astype(str))
        excluded = pd.DataFrame({
            "gene": adata.var_names[~measured],
            "n_samples_unmeasured": (~measurement.measured[np.ix_(present, ~measured)]).sum(axis=0),
            "positioned": positioned[~measured],
        })
        table = write_table(excluded, Path(output_dir) / f"{sample_name}_unmeasured_genes_excluded.tsv")
        gene_metadata = gene_metadata[~gene_metadata["hgnc_symbol"].astype(str).isin(excluded["gene"])]
        plan = replace(plan, measurement={
            "status": "applied", "sample_key": measurement.sample_key, "n_samples": int(len(present)),
            "n_genes_excluded": int(len(excluded)), "n_positioned_genes_excluded": int(excluded["positioned"].sum()),
            "n_positioned_genes_used": int((positioned & measured).sum()), "excluded_genes_table": str(table),
        })
    if not plan.should_run_fastcnv:
        return None, plan

    result = run_fastcnv_pooled_anndata(
        adata,
        gene_metadata,
        sample_key=sample_key,
        layer=layer,
        reference_var=plan.reference_key,
        reference_label=plan.reference_labels,
        config=config or FastCNVConfig(),
        sample_name=sample_name,
        densify_all=h5ad_mode == "dense",
        n_jobs=n_jobs,
        compute_clusters=compute_clusters,
        cell_block_size=cell_block_size,
    )
    write_pooled_outputs(result=result, output_dir=output_dir, sample_name=sample_name, mode="compact")
    attach_tumor_fastcnv_evidence(adata, result, reference_key=plan.reference_key,
                                  all_cells_q75_fallback=all_cells_q75_fallback)
    return result, plan


def attach_tumor_fastcnv_evidence(adata, result: PooledFastCNVResult, *, reference_key: str,
                                  all_cells_q75_fallback: bool = False) -> None:
    meta = result.cell_metadata.reindex(adata.obs_names)
    adata.obs["fastcnv_cnv_fraction"] = meta["cnv_fraction"].astype(float)
    ref_mask = adata.obs[reference_key].astype(str).str.startswith("normal_")
    ref_values = adata.obs.loc[ref_mask, "fastcnv_cnv_fraction"].astype(float)
    if len(ref_values) >= 10:
        med = float(np.nanmedian(ref_values))
        mad = float(np.nanmedian(np.abs(ref_values - med)))
        threshold = max(float(np.nanquantile(ref_values, 0.99)), med + 3.0 * 1.4826 * mad)
        source = "reference_q99_mad"
    elif all_cells_q75_fallback:
        # Opt-in only: a query quantile labels about a quarter of cells tumour-like by construction.
        threshold = float(np.nanquantile(adata.obs["fastcnv_cnv_fraction"].astype(float), 0.75))
        source = "all_cells_q75_fallback"
    else:
        threshold, source = np.nan, "unassessed_insufficient_reference"
    adata.obs["fastcnv_normal_threshold"] = threshold
    adata.obs["fastcnv_threshold_source"] = source
    calls = pd.Series("cnv_indeterminate", index=adata.obs_names, dtype=object)
    calls.loc[ref_mask] = "normal_reference"
    calls.loc[~ref_mask & (adata.obs["fastcnv_cnv_fraction"].astype(float) <= threshold)] = "cnv_normal_like"
    calls.loc[~ref_mask & (adata.obs["fastcnv_cnv_fraction"].astype(float) > threshold)] = "cnv_tumor_like"
    if np.isnan(threshold):
        calls.loc[~ref_mask] = "cnv_unassessed"
    adata.obs["fastcnv_tumor_evidence"] = calls
