from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from fastcnvpy import FastCNVConfig, PooledFastCNVResult, run_fastcnv_pooled_anndata
from fastcnvpy.io import write_pooled_outputs


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
    nonpar = set(nonparenchymal_lineages or DEFAULT_NONPARENCHYMAL_LINEAGES)

    is_parenchymal = major.isin(parenchymal)
    is_normal_nonpar = major.isin(nonpar) | (~is_parenchymal)
    is_normal_par = pd.Series(False, index=obs.index)
    if normal_status_key and normal_status_key in obs.columns:
        is_normal_par = is_parenchymal & obs[normal_status_key].astype(str).isin({str(x) for x in normal_status_labels})

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

    return TumorFastCNVPlan(
        reference_key=reference_key,
        reference_labels=labels,
        n_reference_cells=n_ref,
        n_normal_nonparenchymal=int((ref == "normal_nonparenchymal").sum()),
        n_normal_parenchymal=int((ref == "normal_parenchymal").sum()),
        n_candidate_cells=n_candidate,
        should_run_fastcnv=should_run,
        reason=reason,
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
    h5ad_mode: str = "dense",
    n_jobs: int = 1,
    config: FastCNVConfig | None = None,
    sample_name: str = "tumor_fastcnv",
) -> tuple[PooledFastCNVResult | None, TumorFastCNVPlan]:
    plan = prepare_tumor_fastcnv_reference(
        adata,
        major_lineage_key=major_lineage_key,
        sample_key=sample_key,
        parenchymal_lineages=parenchymal_lineages,
        normal_status_key=normal_status_key,
        normal_status_labels=normal_status_labels,
    )
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
    )
    write_pooled_outputs(result=result, output_dir=output_dir, sample_name=sample_name, mode="compact")
    attach_tumor_fastcnv_evidence(adata, result, reference_key=plan.reference_key)
    return result, plan


def attach_tumor_fastcnv_evidence(adata, result: PooledFastCNVResult, *, reference_key: str) -> None:
    meta = result.cell_metadata.reindex(adata.obs_names)
    adata.obs["fastcnv_cnv_fraction"] = meta["cnv_fraction"].astype(float)
    ref_mask = adata.obs[reference_key].astype(str).str.startswith("normal_")
    ref_values = adata.obs.loc[ref_mask, "fastcnv_cnv_fraction"].astype(float)
    if len(ref_values) >= 10:
        med = float(np.nanmedian(ref_values))
        mad = float(np.nanmedian(np.abs(ref_values - med)))
        threshold = max(float(np.nanquantile(ref_values, 0.99)), med + 3.0 * 1.4826 * mad)
    else:
        threshold = float(np.nanquantile(adata.obs["fastcnv_cnv_fraction"].astype(float), 0.75))
    adata.obs["fastcnv_normal_threshold"] = threshold
    calls = pd.Series("cnv_indeterminate", index=adata.obs_names, dtype=object)
    calls.loc[ref_mask] = "normal_reference"
    calls.loc[~ref_mask & (adata.obs["fastcnv_cnv_fraction"].astype(float) <= threshold)] = "cnv_normal_like"
    calls.loc[~ref_mask & (adata.obs["fastcnv_cnv_fraction"].astype(float) > threshold)] = "cnv_tumor_like"
    adata.obs["fastcnv_tumor_evidence"] = calls
