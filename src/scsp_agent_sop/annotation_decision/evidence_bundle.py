from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from typing import Any, Mapping

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class ClusterEvidence:
    cluster_id: str
    cluster_key: str
    n_cells: int
    n_samples: int | None
    marker_refs: tuple[str, ...]
    top_markers: tuple[str, ...]
    anti_marker_warnings: tuple[str, ...]
    nmf_refs: tuple[str, ...]
    dominant_programs: tuple[str, ...]
    program_warnings: tuple[str, ...]
    reference_refs: tuple[str, ...]
    cnv_refs: tuple[str, ...]
    cnv_summary: dict[str, Any]
    qc_warnings: tuple[str, ...]
    n_donors: int | None
    n_studies: int | None
    nearest_clusters: tuple[str, ...]
    distribution_summary: dict[str, Any]
    program_usage_summary: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return _json_safe(asdict(self))

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ClusterEvidence":
        tuple_fields = {
            "marker_refs",
            "top_markers",
            "anti_marker_warnings",
            "nmf_refs",
            "dominant_programs",
            "program_warnings",
            "reference_refs",
            "cnv_refs",
            "qc_warnings",
            "nearest_clusters",
        }
        normalized = dict(data)
        for field in tuple_fields:
            normalized[field] = tuple(normalized.get(field, ()))
        normalized["cnv_summary"] = dict(normalized.get("cnv_summary", {}))
        normalized["distribution_summary"] = dict(normalized.get("distribution_summary", {}))
        normalized["program_usage_summary"] = dict(normalized.get("program_usage_summary", {}))
        normalized.setdefault("n_donors", None)
        normalized.setdefault("n_studies", None)
        return cls(**normalized)


@dataclass(frozen=True)
class AnnotationEvidenceBundle:
    schema_version: str
    run_id: str
    organism: str
    tissue: str
    disease_context: str | None
    is_tumor: bool
    cluster_key: str
    clusters: tuple[ClusterEvidence, ...]
    file_registry: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["clusters"] = [cluster.to_dict() for cluster in self.clusters]
        return _json_safe(data)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "AnnotationEvidenceBundle":
        normalized = dict(data)
        normalized["clusters"] = tuple(ClusterEvidence.from_dict(item) for item in normalized.get("clusters", ()))
        normalized["file_registry"] = dict(normalized.get("file_registry", {}))
        normalized["is_tumor"] = bool(normalized.get("is_tumor", False))
        return cls(**normalized)

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), ensure_ascii=False, sort_keys=True)

    @classmethod
    def from_json(cls, text: str) -> "AnnotationEvidenceBundle":
        return cls.from_dict(json.loads(text))

    def cluster_by_id(self) -> dict[str, ClusterEvidence]:
        return {cluster.cluster_id: cluster for cluster in self.clusters}


def build_evidence_bundle(
    adata,
    *,
    run_id: str,
    organism: str,
    tissue: str,
    disease_context: str | None,
    is_tumor: bool,
    cluster_key: str,
    file_registry: Mapping[str, Any] | None = None,
    marker_refs: tuple[str, ...] = (),
    reference_refs: tuple[str, ...] = (),
) -> AnnotationEvidenceBundle:
    if cluster_key not in adata.obs:
        raise KeyError(f"cluster_key {cluster_key!r} is not present in adata.obs")
    registry = dict(file_registry or adata.uns.get("file_registry", {}))
    neighbours = _cluster_neighbours(adata, cluster_key)
    clusters = []
    for cluster_id, indices in adata.obs.groupby(cluster_key, observed=True).indices.items():
        obs = adata.obs.iloc[list(indices)]
        clusters.append(
            ClusterEvidence(
                cluster_id=str(cluster_id),
                cluster_key=cluster_key,
                n_cells=int(len(obs)),
                n_samples=_n_unique(obs, "sample_id"),
                marker_refs=tuple(marker_refs),
                top_markers=_extract_top_markers(adata, str(cluster_id)),
                anti_marker_warnings=(),
                nmf_refs=_nmf_refs(adata),
                dominant_programs=_dominant_programs(obs),
                program_warnings=(),
                reference_refs=tuple(reference_refs),
                cnv_refs=_cnv_refs(obs),
                cnv_summary=_cnv_summary(obs),
                qc_warnings=_qc_warnings(obs),
                n_donors=_n_unique(obs, "donor_id"),
                n_studies=_n_unique(obs, "study_id"),
                nearest_clusters=neighbours.get(str(cluster_id), ()),
                distribution_summary=_distribution_summary(obs),
                program_usage_summary=_program_usage_summary(adata, list(indices)),
            )
        )
    return AnnotationEvidenceBundle(
        schema_version="scoop.annotation_evidence_bundle.v1",
        run_id=str(run_id),
        organism=str(organism),
        tissue=str(tissue),
        disease_context=disease_context,
        is_tumor=bool(is_tumor),
        cluster_key=cluster_key,
        clusters=tuple(clusters),
        file_registry=registry,
    )


def _n_unique(obs: pd.DataFrame, key: str) -> int | None:
    if key not in obs:
        return None
    return int(obs[key].nunique())


def _extract_top_markers(adata, cluster_id: str) -> tuple[str, ...]:
    key = f"top_markers_{cluster_id}"
    if key in adata.uns and isinstance(adata.uns[key], (list, tuple, np.ndarray)):
        return tuple(map(str, adata.uns[key]))
    return ()


def _nmf_refs(adata) -> tuple[str, ...]:
    refs = []
    if "X_nmf_usage" in adata.obsm:
        refs.append("obsm.X_nmf_usage")
    if "file_registry" in adata.uns:
        tables = adata.uns["file_registry"].get("tables", {})
        refs.extend(str(key) for key in tables if "program" in str(key).lower() or "nmf" in str(key).lower())
    return tuple(refs)


def _dominant_programs(obs: pd.DataFrame) -> tuple[str, ...]:
    for key in ("dominant_nmf_program", "dominant_program", "program"):
        if key in obs:
            return tuple(map(str, obs[key].dropna().astype(str).value_counts().head(3).index))
    return ()


def _distribution_summary(obs: pd.DataFrame) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for key in ("study_id", "donor_id", "sample_id", "condition", "tissue", "chemistry"):
        if key in obs:
            summary[key] = {
                str(name): int(count) for name, count in obs[key].dropna().astype(str).value_counts().items()
            }
    return summary


def _program_usage_summary(adata, indices: list[int]) -> dict[str, Any]:
    if "X_nmf_usage" not in adata.obsm or not indices:
        return {}
    usage = np.asarray(adata.obsm["X_nmf_usage"])[indices]
    if usage.ndim != 2 or usage.shape[1] == 0:
        return {}
    names = adata.uns.get("nmf_program_names")
    if not isinstance(names, (list, tuple, np.ndarray)) or len(names) != usage.shape[1]:
        names = [f"P{index + 1}" for index in range(usage.shape[1])]
    medians = np.nanmedian(usage, axis=0)
    order = np.argsort(medians)[::-1][: min(5, usage.shape[1])]
    return {
        str(names[index]): {
            "median": float(medians[index]),
            "mean": float(np.nanmean(usage[:, index])),
            "q25": float(np.nanquantile(usage[:, index], 0.25)),
            "q75": float(np.nanquantile(usage[:, index], 0.75)),
        }
        for index in order
    }


def _cluster_neighbours(adata, cluster_key: str, top_n: int = 5) -> dict[str, tuple[str, ...]]:
    graph = None
    for key in ("scvi_connectivities", "connectivities_identity", "connectivities"):
        if key in adata.obsp:
            graph = adata.obsp[key]
            break
    labels = adata.obs[cluster_key].astype(str).to_numpy()
    unique = sorted(set(labels))
    if graph is None:
        return {label: () for label in unique}
    if not hasattr(graph, "tocsr"):
        from scipy import sparse

        graph = sparse.csr_matrix(graph)
    else:
        graph = graph.tocsr()
    result: dict[str, tuple[str, ...]] = {}
    for source in unique:
        source_index = np.flatnonzero(labels == source)
        scores: list[tuple[str, float]] = []
        for target in unique:
            if target == source:
                continue
            target_index = np.flatnonzero(labels == target)
            weight = float(graph[source_index][:, target_index].sum())
            normalizer = max(len(source_index), 1)
            scores.append((target, weight / normalizer))
        scores.sort(key=lambda item: (-item[1], item[0]))
        result[source] = tuple(target for target, weight in scores[:top_n] if weight > 0)
    return result


def _cnv_refs(obs: pd.DataFrame) -> tuple[str, ...]:
    fields = [field for field in _FASTCNV_FIELDS if field in obs]
    return tuple(fields)


_FASTCNV_FIELDS = (
    "fastcnv_reference_pool",
    "fastcnv_cnv_fraction",
    "fastcnv_normal_threshold",
    "fastcnv_tumor_evidence",
)


def _cnv_summary(obs: pd.DataFrame) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    if "fastcnv_tumor_evidence" in obs:
        values = obs["fastcnv_tumor_evidence"].dropna().astype(str)
        summary["fastcnv_tumor_evidence_counts"] = values.value_counts().to_dict()
        summary["has_tumor_evidence"] = bool(
            values.str.lower().isin({"tumor", "malignant", "aneuploid", "positive", "true", "1"}).any()
        )
    if "fastcnv_cnv_fraction" in obs:
        numeric = pd.to_numeric(obs["fastcnv_cnv_fraction"], errors="coerce")
        if numeric.notna().any():
            summary["fastcnv_cnv_fraction_mean"] = float(numeric.mean())
            summary["fastcnv_cnv_fraction_max"] = float(numeric.max())
    if "fastcnv_reference_pool" in obs:
        summary["fastcnv_reference_pool_counts"] = (
            obs["fastcnv_reference_pool"].dropna().astype(str).value_counts().to_dict()
        )
    return _json_safe(summary)


def _qc_warnings(obs: pd.DataFrame) -> tuple[str, ...]:
    warnings = []
    if "qc_pass" in obs and not bool(obs["qc_pass"].astype(bool).all()):
        warnings.append("cluster_contains_qc_failed_cells")
    return tuple(warnings)


def _json_safe(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value
