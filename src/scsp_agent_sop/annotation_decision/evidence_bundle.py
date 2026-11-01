from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field, replace
from typing import Any, Mapping
from pathlib import Path

import numpy as np
import pandas as pd


def bound_evidence_files(bundle, cited_refs) -> list[Path]:
    """Verify referenced registry artifacts before changing annotation outputs."""
    from scsp_agent_sop.storage import sha256_file

    cited_refs = {str(ref).split("::", 1)[0] for ref in cited_refs}
    paths = []
    for category in bundle.file_registry.values():
        if not isinstance(category, Mapping):
            continue
        for key, record in category.items():
            if key in cited_refs and isinstance(record, Mapping) and record.get("path"):
                path = Path(record["path"]).resolve()
                expected = record.get("sha256") or bundle.deg_tables.get(key, {}).get("sha256")
                if not path.is_file() or (expected and sha256_file(path) != expected):
                    raise ValueError(f"Cited annotation evidence is missing or changed: {path}")
                paths.append(path)
    return paths


def partition_sha256(adata, cluster_key: str) -> str:
    """Bind the decision units to the exact obs axis and unit membership."""
    digest = hashlib.sha256(f"{cluster_key}\0".encode())
    digest.update("\0".join(adata.obs_names.astype(str)).encode())
    digest.update(b"\1")
    digest.update("\0".join(adata.obs[cluster_key].astype(str)).encode())
    return digest.hexdigest()


DEG_GROUP_COLUMNS = ("group", "cluster", "cell_type")
DEG_GENE_COLUMNS = ("names", "gene", "feature")
DEG_SCORE_COLUMNS = ("scores", "score", "z")
TOP_MARKERS_FROM_TABLE = 25


def deg_tables_from_registry(registry: Mapping[str, Any]) -> tuple[dict[str, dict[str, Any]], dict[str, pd.DataFrame]]:
    """Read registered DEG tables: which partition (label_key) and groups they test.

    A table is a candidate when its record declares label_key or its schema
    names markers; it is a DEG table only if it has group and gene columns.
    Unreadable files are left out, so they can never support a citation.
    """
    from scsp_agent_sop.storage import sha256_file

    tables, frames = {}, {}
    for key, record in dict(registry.get("tables", {})).items():
        if not isinstance(record, Mapping) or not record.get("path"):
            continue
        if "label_key" not in record and "markers" not in str(record.get("schema", "")):
            continue
        path = Path(record["path"])
        if not path.is_file():
            continue
        if path.suffix.lower() == ".parquet":
            import pyarrow.parquet as pq

            names = pq.read_schema(path).names
            wanted = (*DEG_GROUP_COLUMNS, *DEG_GENE_COLUMNS, *DEG_SCORE_COLUMNS, "label_key")
            frame = pd.read_parquet(path, columns=[name for name in wanted if name in names])
        else:
            # Groups are labels ("01" is not 1); scores are numbers.
            frame = pd.read_csv(path, sep="," if path.suffix.lower() == ".csv" else "\t", dtype=str)
            for name in DEG_SCORE_COLUMNS:
                if name in frame:
                    frame[name] = pd.to_numeric(frame[name], errors="coerce")
        group = next((name for name in DEG_GROUP_COLUMNS if name in frame), None)
        gene = next((name for name in DEG_GENE_COLUMNS if name in frame), None)
        if group is None or gene is None:
            continue
        label_key = record.get("label_key")
        if label_key is None and "label_key" in frame and frame["label_key"].nunique() == 1:
            label_key = str(frame["label_key"].iloc[0])
        frame = frame.assign(**{group: frame[group].astype(str)})
        tables[str(key)] = {"label_key": None if label_key is None else str(label_key),
                            "groups": sorted(frame[group].unique().tolist()),
                            "sha256": record.get("sha256") or sha256_file(path)}
        score = next((name for name in DEG_SCORE_COLUMNS if name in frame), None)
        if score is not None:
            frame = frame.sort_values(score, ascending=False, kind="stable")
        frames[str(key)] = frame.rename(columns={group: "_group", gene: "_gene"})[["_group", "_gene"]]
    return tables, frames


def bind_legacy_bundle(bundle: "AnnotationEvidenceBundle", adata) -> "AnnotationEvidenceBundle":
    """Rebuild content bindings a v1 bundle lacks; the commit reports this as limited."""
    deg_tables, _ = deg_tables_from_registry(bundle.file_registry)
    return replace(bundle, partition_sha256=partition_sha256(adata, bundle.cluster_key), deg_tables=deg_tables)


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
    # n_samples counts samples with >= min_cells_per_sample cells; this is the raw count.
    n_samples_raw: int | None = None
    # DEG tables testing this partition (label_key == cluster_key) that contain this group.
    deg_refs: tuple[str, ...] = ()

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
            "deg_refs",
        }
        normalized = dict(data)
        for name in tuple_fields:
            normalized[name] = tuple(normalized.get(name, ()))
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
    # Absent in v1 bundles; bind_legacy_bundle rebuilds them from the input.
    partition_sha256: str | None = None
    deg_tables: dict[str, dict[str, Any]] = field(default_factory=dict)
    composition_keys: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["clusters"] = [cluster.to_dict() for cluster in self.clusters]
        return _json_safe(data)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "AnnotationEvidenceBundle":
        normalized = dict(data)
        normalized["clusters"] = tuple(ClusterEvidence.from_dict(item) for item in normalized.get("clusters", ()))
        normalized["file_registry"] = dict(normalized.get("file_registry", {}))
        normalized["deg_tables"] = dict(normalized.get("deg_tables") or {})
        normalized["composition_keys"] = dict(normalized.get("composition_keys") or {})
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
    sample_key: str = "sample_id",
    donor_key: str = "donor_id",
    study_key: str = "study_id",
    min_cells_per_sample: int = 20,
) -> AnnotationEvidenceBundle:
    if cluster_key not in adata.obs:
        raise KeyError(f"cluster_key {cluster_key!r} is not present in adata.obs")
    registry = dict(file_registry or adata.uns.get("file_registry", {}))
    if not marker_refs:
        marker_refs = tuple(str(key) for key, record in registry.get("tables", {}).items()
                            if isinstance(record, Mapping) and str(record.get("schema", "")).startswith("cluster_markers"))
    deg_tables, deg_frames = deg_tables_from_registry(registry)
    keyed_deg = any(table["label_key"] is not None for table in deg_tables.values())
    composition_keys = {"sample_key": sample_key, "donor_key": donor_key, "study_key": study_key,
                        "min_cells_per_sample": int(min_cells_per_sample)}
    neighbours = _cluster_neighbours(adata, cluster_key)
    clusters = []
    for cluster_id, indices in adata.obs.groupby(cluster_key, observed=True).indices.items():
        obs = adata.obs.iloc[list(indices)]
        cluster_id = str(cluster_id)
        deg_refs = tuple(ref for ref, table in deg_tables.items()
                         if table["label_key"] == cluster_key and cluster_id in table["groups"])
        # A registered marker table that tests another partition, or lacks this
        # group, must not look like this unit's DEG (FC-06, C2-02).
        unit_marker_refs = tuple(ref for ref in marker_refs if ref not in deg_tables or ref in deg_refs)
        # uns top_markers_<id> carry no partition; they are a fallback only for
        # registries that predate label-keyed DEG tables.
        top_markers = () if keyed_deg else _extract_top_markers(adata, cluster_id)
        if deg_refs:
            frame = deg_frames[deg_refs[0]]
            top_markers = tuple(frame.loc[frame["_group"] == cluster_id, "_gene"].astype(str).head(TOP_MARKERS_FROM_TABLE))
        clusters.append(
            ClusterEvidence(
                cluster_id=cluster_id,
                cluster_key=cluster_key,
                n_cells=int(len(obs)),
                n_samples=_n_unique(obs, sample_key, min_cells_per_sample),
                marker_refs=unit_marker_refs,
                top_markers=top_markers,
                anti_marker_warnings=(),
                nmf_refs=_nmf_refs(adata),
                dominant_programs=_dominant_programs(obs),
                program_warnings=(),
                reference_refs=tuple(reference_refs),
                cnv_refs=_cnv_refs(obs),
                cnv_summary=_cnv_summary(obs),
                qc_warnings=_qc_warnings(obs),
                n_donors=_n_unique(obs, donor_key),
                n_studies=_n_unique(obs, study_key),
                nearest_clusters=neighbours.get(cluster_id, ()),
                distribution_summary=_distribution_summary(obs, (study_key, donor_key, sample_key)),
                program_usage_summary=_program_usage_summary(adata, list(indices)),
                n_samples_raw=_n_unique(obs, sample_key),
                deg_refs=deg_refs,
            )
        )
    return AnnotationEvidenceBundle(
        schema_version="scoop.annotation_evidence_bundle.v2",
        run_id=str(run_id),
        organism=str(organism),
        tissue=str(tissue),
        disease_context=disease_context,
        is_tumor=bool(is_tumor),
        cluster_key=cluster_key,
        clusters=tuple(clusters),
        file_registry=registry,
        partition_sha256=partition_sha256(adata, cluster_key),
        deg_tables=deg_tables,
        composition_keys=composition_keys,
    )


def _n_unique(obs: pd.DataFrame, key: str, min_cells: int = 1) -> int | None:
    if key not in obs:
        return None
    return int((obs[key].value_counts() >= max(min_cells, 1)).sum())


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


def _distribution_summary(obs: pd.DataFrame, composition_keys: tuple[str, ...]) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for key in dict.fromkeys((*composition_keys, "condition", "tissue", "chemistry")):
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
    "fastcnv_threshold_source",
)
# Calls without a reference threshold: not assessed, never CNV-negative or normal.
FASTCNV_UNASSESSED_CALLS = {"cnv_unassessed"}


def _cnv_summary(obs: pd.DataFrame) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    if "fastcnv_tumor_evidence" in obs:
        values = obs["fastcnv_tumor_evidence"].dropna().astype(str)
        summary["fastcnv_tumor_evidence_counts"] = values.value_counts().to_dict()
        assessed = values[~values.isin(FASTCNV_UNASSESSED_CALLS)]
        summary["fastcnv_assessed_cells"] = int(len(assessed))
        summary["fastcnv_unassessed_cells"] = int(len(values) - len(assessed))
        # None, not False, when no cell was assessed.
        summary["has_tumor_evidence"] = bool(
            assessed.str.lower().isin({"cnv_tumor_like", "tumor", "malignant", "aneuploid", "positive", "true", "1"}).any()
        ) if len(assessed) else None
    if "fastcnv_threshold_source" in obs:
        summary["fastcnv_threshold_source_counts"] = (
            obs["fastcnv_threshold_source"].dropna().astype(str).value_counts().to_dict()
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
