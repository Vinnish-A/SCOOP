from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping
import json

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

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

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
        }
        normalized = dict(data)
        for field in tuple_fields:
            normalized[field] = tuple(normalized.get(field, ()))
        normalized["cnv_summary"] = dict(normalized.get("cnv_summary", {}))
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
        return data

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
    if key in adata.uns and isinstance(adata.uns[key], (list, tuple)):
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


def _cnv_refs(obs: pd.DataFrame) -> tuple[str, ...]:
    fields = [field for field in _FASTCNV_FIELDS if field in obs]
    return tuple(fields)


_FASTCNV_FIELDS = ("fastcnv_reference_pool", "fastcnv_cnv_fraction", "fastcnv_normal_threshold", "fastcnv_tumor_evidence")


def _cnv_summary(obs: pd.DataFrame) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    if "fastcnv_tumor_evidence" in obs:
        values = obs["fastcnv_tumor_evidence"].dropna().astype(str)
        summary["fastcnv_tumor_evidence_counts"] = values.value_counts().to_dict()
        summary["has_tumor_evidence"] = bool(values.str.lower().isin({"tumor", "malignant", "aneuploid", "positive", "true", "1"}).any())
    if "fastcnv_cnv_fraction" in obs:
        numeric = pd.to_numeric(obs["fastcnv_cnv_fraction"], errors="coerce")
        if numeric.notna().any():
            summary["fastcnv_cnv_fraction_mean"] = float(numeric.mean())
            summary["fastcnv_cnv_fraction_max"] = float(numeric.max())
    if "fastcnv_reference_pool" in obs:
        summary["fastcnv_reference_pool_counts"] = obs["fastcnv_reference_pool"].dropna().astype(str).value_counts().to_dict()
    return _json_safe(summary)


def _qc_warnings(obs: pd.DataFrame) -> tuple[str, ...]:
    warnings = []
    if "qc_pass" in obs and not bool(obs["qc_pass"].astype(bool).all()):
        warnings.append("cluster_contains_qc_failed_cells")
    return tuple(warnings)


def _json_safe(value: Any) -> Any:
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
