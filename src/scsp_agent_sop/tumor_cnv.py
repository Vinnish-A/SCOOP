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
    # "default" (DEFAULT_NONPARENCHYMAL_LINEAGES), "explicit", or
    # "not_applied_reference_eligible_key" when a review mask selects references (DC4).
    lineage_vocabulary: str = "default"
    min_reference_cells: int = 20
    min_reference_cells_per_sample: int = 5
    # Reference cells per sample that the engine uses for that sample's own
    # baseline; the other samples are scored against the pooled baseline (DV7).
    reference_cells_by_sample: dict[str, int] = field(default_factory=dict)
    samples_without_reference: list[str] = field(default_factory=list)
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
    min_reference_cells_per_sample: int = 5,
    all_cells_q75_fallback: bool = False,
) -> TumorFastCNVPlan:
    """Build the pooled-reference labels used by the tumor FastCNV step.

    ``min_reference_cells`` pooled references give a reference-derived
    threshold; with fewer, the plan runs only when ``all_cells_q75_fallback``
    is opted in. A sample keeps its own baseline when one pool label has
    ``min_reference_cells_per_sample`` cells in it (the engine's rule); at least
    one sample must. With ``reference_eligible_key`` and no explicit
    ``nonparenchymal_lineages``, every eligible non-parenchymal cell is a
    reference: the review mask, not a built-in vocabulary, decides (DC4).
    """

    if major_lineage_key not in adata.obs.columns:
        raise ValueError(f"major_lineage_key {major_lineage_key!r} is not present in obs")
    if sample_key not in adata.obs.columns:
        raise ValueError(f"sample_key {sample_key!r} is not present in obs")

    obs = adata.obs
    major = obs[major_lineage_key].astype(str)
    parenchymal = set(map(str, parenchymal_lineages))
    if nonparenchymal_lineages is not None:
        nonpar, vocabulary = set(nonparenchymal_lineages), "explicit"
    elif reference_eligible_key:
        nonpar, vocabulary = None, "not_applied_reference_eligible_key"
    else:
        nonpar, vocabulary = set(DEFAULT_NONPARENCHYMAL_LINEAGES), "default"
    in_vocabulary = major.isin(nonpar) if nonpar is not None else pd.Series(True, index=obs.index)
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
    is_normal_nonpar = in_vocabulary & ~is_parenchymal & eligible
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
    is_ref = ref.isin(labels)
    n_ref = int(is_ref.sum())
    n_candidate = int((ref == "candidate_tumor_or_mixed").sum())
    per_label = pd.crosstab(obs.loc[is_ref, sample_key], ref[is_ref])  # samples x labels; NaN samples drop out
    own = per_label.where(per_label >= min_reference_cells_per_sample, 0).sum(axis=1)
    reference_cells_by_sample = {str(sample): int(n) for sample, n in own.items() if n > 0}
    samples = sorted(set(obs[sample_key].dropna().astype(str)))
    enough_pooled = n_ref >= min_reference_cells
    should_run = bool(n_candidate > 0 and reference_cells_by_sample and (enough_pooled or all_cells_q75_fallback))
    if n_candidate == 0:
        reason = "no candidate parenchymal/tumor cells"
    elif not reference_cells_by_sample:
        reason = (f"no sample has {min_reference_cells_per_sample} reference cells of one pool label "
                  "(tumor_fastcnv.min_reference_cells_per_sample)")
    elif not should_run:
        reason = (f"{n_ref} pooled reference cells, below tumor_fastcnv.min_reference_cells={min_reference_cells}; "
                  "opt in to all_cells_q75_fallback for a query-quantile threshold")
    elif int(is_normal_par.sum()) > 0:
        reason = "normal parenchymal and tumor/candidate cells are separable; FastCNV will provide orthogonal CNV evidence"
    else:
        reason = "normal parenchymal is mixed or unavailable; FastCNV will identify CNV-supported tumor cells using pooled nonparenchymal reference"
    if should_run and not enough_pooled:
        reason += f"; only {n_ref} pooled reference cells, so the threshold is the opted-in all-cells q75"
    # Reviewed nonparenchymal cells outside the vocabulary are excluded, not hidden.
    dropped = major[eligible & ~is_parenchymal & ~in_vocabulary].value_counts()

    return TumorFastCNVPlan(
        reference_key=reference_key,
        reference_labels=labels,
        n_reference_cells=n_ref,
        n_normal_nonparenchymal=int((ref == "normal_nonparenchymal").sum()),
        n_normal_parenchymal=int((ref == "normal_parenchymal").sum()),
        n_candidate_cells=n_candidate,
        should_run_fastcnv=should_run,
        reason=reason,
        nonparenchymal_lineages=sorted(nonpar) if nonpar is not None else [],
        dropped_reference_by_lineage={str(k): int(v) for k, v in dropped.items()},
        lineage_vocabulary=vocabulary,
        min_reference_cells=int(min_reference_cells),
        min_reference_cells_per_sample=int(min_reference_cells_per_sample),
        reference_cells_by_sample=reference_cells_by_sample,
        samples_without_reference=[sample for sample in samples if sample not in reference_cells_by_sample],
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
    min_reference_cells: int = 20,
    min_reference_cells_per_sample: int = 5,
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
    ``plan.measurement``. ``obs["fastcnv_sample_has_reference"]`` is False
    for cells whose sample had no own baseline (``plan.samples_without_reference``).
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
        min_reference_cells=min_reference_cells,
        min_reference_cells_per_sample=min_reference_cells_per_sample,
        all_cells_q75_fallback=all_cells_q75_fallback,
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
        min_reference_cells_per_sample=plan.min_reference_cells_per_sample,
    )
    write_pooled_outputs(result=result, output_dir=output_dir, sample_name=sample_name, mode="compact")
    attach_tumor_fastcnv_evidence(adata, result, reference_key=plan.reference_key,
                                  all_cells_q75_fallback=all_cells_q75_fallback, min_reference_cells=min_reference_cells)
    samples = adata.obs[sample_key]
    adata.obs["fastcnv_sample_has_reference"] = (
        samples.notna() & samples.astype(str).isin(list(plan.reference_cells_by_sample))).to_numpy()
    return result, plan


def attach_tumor_fastcnv_evidence(adata, result: PooledFastCNVResult, *, reference_key: str,
                                  all_cells_q75_fallback: bool = False, min_reference_cells: int = 20) -> None:
    meta = result.cell_metadata.reindex(adata.obs_names)
    adata.obs["fastcnv_cnv_fraction"] = meta["cnv_fraction"].astype(float)
    ref_mask = adata.obs[reference_key].astype(str).str.startswith("normal_")
    # A reference cell the engine did not score (NaN) is no evidence for the threshold.
    ref_values = adata.obs.loc[ref_mask, "fastcnv_cnv_fraction"].astype(float).dropna()
    if len(ref_values) >= min_reference_cells:
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
