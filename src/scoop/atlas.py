"""Executable large-atlas integration, annotation reporting, and scPoli gate."""

from __future__ import annotations

import importlib.metadata
import json
from datetime import datetime, timezone
from html import escape
from pathlib import Path
from typing import Any, Iterable

import anndata as ad
import numpy as np
import pandas as pd

from scsp_agent_sop.annotation import run_markers_scanpy
from scsp_agent_sop.annotation_decision.evidence_bundle import build_evidence_bundle
from scsp_agent_sop.config import deep_get, read_yaml, resolve_run_root
from scsp_agent_sop.storage import init_file_registry, register_file, write_table

from .approval import verify as verify_approval
from .io import atomic_json, atomic_text

ATLAS_REPORT_SCHEMA = "scoop.atlas-report.v1"
LABEL_COLUMNS = (
    "cluster_id",
    "cell_type_lvl1",
    "cell_type_lvl2",
    "cell_type_lvl3",
    "cell_state",
    "scvi_curated_cell_type",
    "annotation_confidence",
    "annotation_status",
    "annotation_source",
    "annotation_version",
    "annotation_rationale",
)


def _version(distribution: str) -> str:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return "unavailable"


def _require_obs(adata: ad.AnnData, keys: Iterable[str]) -> None:
    missing = [key for key in keys if key not in adata.obs]
    if missing:
        raise KeyError(f"required atlas obs fields are missing: {', '.join(missing)}")


def integrate(
    *,
    config_path: str | Path,
    input_path: str | Path | None = None,
    output_path: str | Path | None = None,
    max_epochs: int | None = None,
) -> Path:
    try:
        import scanpy as sc
        import scvi
    except ImportError as exc:
        raise RuntimeError("atlas integration requires scvi-tools in the managed SCOOP runtime") from exc

    config_path = Path(config_path).resolve()
    config = read_yaml(config_path)
    run_root = resolve_run_root(config_path, config).resolve()
    source = Path(input_path) if input_path else run_root / "artifacts" / "adata_qc.h5ad"
    target = Path(output_path) if output_path else run_root / "artifacts" / "adata_scvi.h5ad"
    adata = ad.read_h5ad(source)
    init_file_registry(adata, deep_get(config, "run.run_id", run_root.name))

    counts_layer = str(deep_get(config, "atlas.scvi.counts_layer", "counts"))
    batch_key = str(deep_get(config, "atlas.scvi.batch_key", deep_get(config, "keys.sample", "sample_id")))
    requested_categorical = list(deep_get(config, "atlas.scvi.categorical_covariate_keys", []))
    categorical = [key for key in requested_categorical if key in adata.obs]
    _require_obs(adata, [batch_key])
    if counts_layer not in adata.layers:
        raise KeyError(f"scVI requires raw counts in adata.layers[{counts_layer!r}]")

    scvi.settings.seed = int(deep_get(config, "run.random_seed", 0))
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=counts_layer,
        batch_key=batch_key,
        categorical_covariate_keys=categorical or None,
    )
    model = scvi.model.SCVI(
        adata,
        n_latent=int(deep_get(config, "atlas.scvi.n_latent", 30)),
        gene_likelihood=str(deep_get(config, "atlas.scvi.gene_likelihood", "nb")),
    )
    epochs = int(max_epochs or deep_get(config, "atlas.scvi.max_epochs", 400))
    model.train(max_epochs=epochs, early_stopping=True)
    adata.obsm["X_scVI"] = model.get_latent_representation()
    neighbours = int(deep_get(config, "atlas.graph.n_neighbors", 30))
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=neighbours, key_added="scvi")
    cluster_key = str(deep_get(config, "atlas.graph.cluster_key", "leiden_scvi"))
    sc.tl.leiden(
        adata,
        neighbors_key="scvi",
        key_added=cluster_key,
        resolution=float(deep_get(config, "atlas.graph.leiden_resolution", 1.0)),
        random_state=int(deep_get(config, "run.random_seed", 0)),
    )
    sc.tl.umap(adata, neighbors_key="scvi", random_state=int(deep_get(config, "run.random_seed", 0)))

    model_path = run_root / "atlas" / "models" / "scvi"
    if model_path.exists() and any(model_path.iterdir()):
        raise FileExistsError(f"refusing to overwrite existing scVI model: {model_path}")
    model_path.parent.mkdir(parents=True, exist_ok=True)
    model.save(model_path, overwrite=False, save_anndata=False)
    provenance = {
        "schema_version": "scoop.scvi-run.v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "input": str(source.resolve()),
        "output": str(target.resolve()),
        "model": str(model_path.resolve()),
        "versions": {"scvi-tools": _version("scvi-tools"), "torch": _version("torch")},
        "parameters": {
            "counts_layer": counts_layer,
            "batch_key": batch_key,
            "categorical_covariate_keys": categorical,
            "skipped_categorical_covariate_keys": sorted(set(requested_categorical) - set(categorical)),
            "n_latent": int(deep_get(config, "atlas.scvi.n_latent", 30)),
            "max_epochs": epochs,
            "n_neighbors": neighbours,
            "cluster_key": cluster_key,
            "leiden_resolution": float(deep_get(config, "atlas.graph.leiden_resolution", 1.0)),
        },
        "annotation_contract": (
            "labels are produced by the downstream annotation step in X_scVI, not by the scVI model head"
        ),
    }
    provenance_path = atomic_json(run_root / "atlas" / "scvi-run.json", provenance)
    register_file(adata, key="scvi_model", path=model_path, category="artifacts", schema="scvi.model")
    register_file(adata, key="scvi_run", path=provenance_path, category="reports", schema="scoop.scvi-run.v1")
    target.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(target)
    return target


def _read_labels(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, sep="\t" if path.suffix.lower() in {".tsv", ".txt"} else ",")


def _label_template(adata: ad.AnnData, cluster_key: str, version: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for cluster_id in sorted(adata.obs[cluster_key].astype(str).unique()):
        subset = adata.obs.loc[adata.obs[cluster_key].astype(str).eq(cluster_id)]
        row: dict[str, Any] = {
            "cluster_id": cluster_id,
            "cell_type_lvl1": "Unresolved",
            "cell_type_lvl2": "Unresolved",
            "cell_type_lvl3": "Unresolved",
            "cell_state": "",
            "scvi_curated_cell_type": f"Unresolved_{cluster_id}",
            "annotation_confidence": "low",
            "annotation_status": "unlabeled",
            "annotation_source": "scvi_cluster_curated",
            "annotation_version": version,
            "annotation_rationale": "Requires marker, reference, neighbour, distribution, and programme review.",
        }
        for key in LABEL_COLUMNS[1:6]:
            if key in subset and subset[key].notna().any() and subset[key].astype(str).nunique() == 1:
                row[key] = str(subset[key].dropna().iloc[0])
        rows.append(row)
    return pd.DataFrame(rows, columns=LABEL_COLUMNS)


def annotate(
    *,
    config_path: str | Path,
    input_path: str | Path | None = None,
    output_path: str | Path | None = None,
    labels_path: str | Path | None = None,
) -> tuple[Path, Path, Path]:
    config_path = Path(config_path).resolve()
    config = read_yaml(config_path)
    run_root = resolve_run_root(config_path, config).resolve()
    source = Path(input_path) if input_path else run_root / "artifacts" / "adata_scvi.h5ad"
    target = Path(output_path) if output_path else run_root / "artifacts" / "adata_scvi_annotated.h5ad"
    adata = ad.read_h5ad(source)
    init_file_registry(adata, deep_get(config, "run.run_id", run_root.name))
    cluster_key = str(deep_get(config, "atlas.graph.cluster_key", "leiden_scvi"))
    _require_obs(adata, [cluster_key])

    marker_layer = str(deep_get(config, "atlas.annotation.marker_layer", "log1p_norm"))
    if marker_layer not in adata.layers:
        marker_layer = "__adata_X__"
    markers = run_markers_scanpy(
        adata,
        groupby=cluster_key,
        layer=marker_layer,
        method=str(deep_get(config, "annotation.default_marker_method", "wilcoxon")),
        n_genes=int(deep_get(config, "annotation.marker_n_genes", 100)),
    )
    marker_path = write_table(markers, run_root / "atlas" / "tables" / "scvi_cluster_markers.parquet")
    group_column = next((key for key in ("group", "cluster", "cell_type") if key in markers), None)
    gene_column = next((key for key in ("names", "gene", "feature") if key in markers), None)
    if group_column and gene_column:
        for cluster_id, frame in markers.groupby(group_column, observed=True):
            adata.uns[f"top_markers_{cluster_id}"] = frame[gene_column].dropna().astype(str).head(100).tolist()

    version = str(deep_get(config, "atlas.annotation.version", "atlas-v1"))
    labels = _read_labels(Path(labels_path)) if labels_path else _label_template(adata, cluster_key, version)
    missing = [key for key in LABEL_COLUMNS if key not in labels]
    if missing:
        raise ValueError(f"label table is missing required columns: {', '.join(missing)}")
    duplicate = labels["cluster_id"].astype(str).duplicated()
    if duplicate.any():
        raise ValueError("label table contains duplicate cluster_id values")
    labels["cluster_id"] = labels["cluster_id"].astype(str)
    cluster_ids = set(adata.obs[cluster_key].astype(str))
    if cluster_ids != set(labels["cluster_id"]):
        raise ValueError("label table cluster_id values must exactly match the scVI clusters")
    index = labels.set_index("cluster_id")
    cluster_values = adata.obs[cluster_key].astype(str)
    for key in LABEL_COLUMNS[1:]:
        adata.obs[key] = cluster_values.map(index[key]).astype(str).to_numpy()

    marker_refs = tuple(str(item) for item in deep_get(config, "atlas.annotation.marker_references", []))
    reference_refs = tuple(str(item) for item in deep_get(config, "atlas.annotation.external_references", []))
    bundle = build_evidence_bundle(
        adata,
        run_id=str(deep_get(config, "run.run_id", run_root.name)),
        organism=str(deep_get(config, "run.organism", "unknown")),
        tissue=str(deep_get(config, "run.tissue", "unknown")),
        disease_context=deep_get(config, "run.disease_context"),
        is_tumor=bool(deep_get(config, "run.is_tumor", False)),
        cluster_key=cluster_key,
        marker_refs=marker_refs,
        reference_refs=reference_refs,
    )
    label_path = write_table(labels, run_root / "atlas" / "tables" / "scvi_curated_labels.tsv")
    evidence_path = atomic_json(run_root / "atlas" / "tables" / "scvi_annotation_evidence.json", bundle.to_dict())
    register_file(adata, key="scvi_cluster_markers", path=marker_path, schema="cluster_markers.v1")
    register_file(adata, key="scvi_curated_labels", path=label_path, schema="scoop.atlas-labels.v1")
    register_file(
        adata, key="scvi_annotation_evidence", path=evidence_path, schema="scoop.annotation_evidence_bundle.v1"
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(target)
    return target, label_path, evidence_path


def _distribution_table(adata: ad.AnnData, cluster_key: str) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for key in ("study_id", "donor_id", "sample_id", "condition", "tissue", "chemistry"):
        if key not in adata.obs:
            continue
        table = (
            adata.obs.assign(_cluster=adata.obs[cluster_key].astype(str), _group=adata.obs[key].astype(str))
            .groupby(["_cluster", "_group"], observed=True)
            .size()
            .rename("n_cells")
            .reset_index()
            .rename(columns={"_cluster": "cluster_id", "_group": "group"})
        )
        table["metadata_key"] = key
        table["p_group_given_subtype"] = table["n_cells"] / table.groupby("cluster_id")["n_cells"].transform("sum")
        table["p_subtype_given_group"] = table["n_cells"] / table.groupby("group")["n_cells"].transform("sum")
        frames.append(table)
    columns = ["cluster_id", "metadata_key", "group", "n_cells", "p_group_given_subtype", "p_subtype_given_group"]
    return pd.concat(frames, ignore_index=True)[columns] if frames else pd.DataFrame(columns=columns)


def _neighbour_table(adata: ad.AnnData, cluster_key: str, top_n: int = 5) -> pd.DataFrame:
    graph = next(
        (
            adata.obsp[key]
            for key in ("scvi_connectivities", "connectivities_identity", "connectivities")
            if key in adata.obsp
        ),
        None,
    )
    labels = adata.obs[cluster_key].astype(str).to_numpy()
    unique = sorted(set(labels))
    latent = np.asarray(adata.obsm["X_scVI"]) if "X_scVI" in adata.obsm else None
    rows: list[dict[str, Any]] = []
    if graph is not None:
        from scipy import sparse

        graph = sparse.csr_matrix(graph)
    for source in unique:
        source_index = np.flatnonzero(labels == source)
        candidates = []
        for target in unique:
            if target == source:
                continue
            target_index = np.flatnonzero(labels == target)
            edge_weight = (
                float(graph[source_index][:, target_index].sum() / max(len(source_index), 1))
                if graph is not None
                else 0.0
            )
            centroid_distance = None
            if latent is not None:
                centroid_distance = float(
                    np.linalg.norm(latent[source_index].mean(axis=0) - latent[target_index].mean(axis=0))
                )
            candidates.append((target, edge_weight, centroid_distance))
        candidates.sort(key=lambda item: (-item[1], item[2] if item[2] is not None else float("inf"), item[0]))
        for rank, (target, edge_weight, distance) in enumerate(candidates[:top_n], start=1):
            rows.append(
                {
                    "source_cluster": source,
                    "target_cluster": target,
                    "rank": rank,
                    "directed_knn_edge_weight": edge_weight,
                    "latent_centroid_distance": distance,
                    "review_decision": "uncertain",
                    "review_rationale": "Compare distinguishing markers, composition, and cNMF usage.",
                }
            )
    return pd.DataFrame(
        rows,
        columns=[
            "source_cluster",
            "target_cluster",
            "rank",
            "directed_knn_edge_weight",
            "latent_centroid_distance",
            "review_decision",
            "review_rationale",
        ],
    )


def _programme_table(bundle) -> pd.DataFrame:
    rows = []
    for cluster in bundle.clusters:
        for name, values in cluster.program_usage_summary.items():
            rows.append({"cluster_id": cluster.cluster_id, "programme": name, **values})
    return pd.DataFrame(rows, columns=["cluster_id", "programme", "median", "mean", "q25", "q75"])


def _programme_classification(programmes: pd.DataFrame) -> pd.DataFrame:
    names = sorted(programmes["programme"].astype(str).unique()) if not programmes.empty else []
    nuisance_terms = ("stress", "ribo", "mito", "cycle", "dissociation", "ambient", "hypoxia")
    rows = []
    for name in names:
        nuisance = any(term in name.lower() for term in nuisance_terms)
        rows.append(
            {
                "programme": name,
                "programme_class": "technical_or_nuisance" if nuisance else "mixed_or_unresolved",
                "classification_rationale": "Name matches a nuisance/state programme"
                if nuisance
                else "Requires top-gene biological review",
            }
        )
    return pd.DataFrame(rows, columns=["programme", "programme_class", "classification_rationale"])


def _reference_concordance(adata: ad.AnnData, cluster_key: str) -> pd.DataFrame:
    reference_keys = [
        key
        for key in adata.obs
        if key.startswith("reference_") and (key.endswith("label") or key.endswith("prediction"))
    ]
    rows = []
    for key in reference_keys:
        for cluster, frame in adata.obs.groupby(cluster_key, observed=True):
            values = frame[key].dropna().astype(str)
            if values.empty:
                continue
            counts = values.value_counts()
            rows.append(
                {
                    "cluster_id": str(cluster),
                    "reference_key": key,
                    "top_reference_label": str(counts.index[0]),
                    "top_label_fraction": float(counts.iloc[0] / counts.sum()),
                    "n_reference_labeled": int(counts.sum()),
                }
            )
    return pd.DataFrame(
        rows,
        columns=["cluster_id", "reference_key", "top_reference_label", "top_label_fraction", "n_reference_labeled"],
    )


def report(
    *,
    config_path: str | Path,
    input_path: str | Path | None = None,
    output_path: str | Path | None = None,
) -> tuple[Path, Path, Path]:
    config_path = Path(config_path).resolve()
    config = read_yaml(config_path)
    run_root = resolve_run_root(config_path, config).resolve()
    source = Path(input_path) if input_path else run_root / "artifacts" / "adata_scvi_annotated.h5ad"
    target = Path(output_path) if output_path else run_root / "reports" / "scvi_annotation_report.html"
    adata = ad.read_h5ad(source)
    cluster_key = str(deep_get(config, "atlas.graph.cluster_key", "leiden_scvi"))
    _require_obs(adata, [cluster_key, "scvi_curated_cell_type", "annotation_confidence", "annotation_status"])
    bundle = build_evidence_bundle(
        adata,
        run_id=str(deep_get(config, "run.run_id", run_root.name)),
        organism=str(deep_get(config, "run.organism", "unknown")),
        tissue=str(deep_get(config, "run.tissue", "unknown")),
        disease_context=deep_get(config, "run.disease_context"),
        is_tumor=bool(deep_get(config, "run.is_tumor", False)),
        cluster_key=cluster_key,
        marker_refs=tuple(str(item) for item in deep_get(config, "atlas.annotation.marker_references", [])),
        reference_refs=tuple(str(item) for item in deep_get(config, "atlas.annotation.external_references", [])),
    )
    table_dir = run_root / "atlas" / "tables"
    distributions = _distribution_table(adata, cluster_key)
    neighbours = _neighbour_table(adata, cluster_key)
    programmes = _programme_table(bundle)
    programme_classes = _programme_classification(programmes)
    reference_concordance = _reference_concordance(adata, cluster_key)
    distribution_path = write_table(distributions, table_dir / "scvi_subgroup_distribution.parquet")
    neighbour_path = write_table(neighbours, table_dir / "scvi_cluster_neighbours.tsv")
    programme_path = write_table(programmes, table_dir / "cnmf_usage_by_subgroup.parquet")
    programme_class_path = write_table(programme_classes, table_dir / "cnmf_program_classification.tsv")
    reference_path = write_table(reference_concordance, table_dir / "reference_mapping_concordance.tsv")
    label_path = table_dir / "scvi_curated_labels.tsv"
    if not label_path.is_file():
        labels = _label_template(adata, cluster_key, str(deep_get(config, "atlas.annotation.version", "atlas-v1")))
        label_path = write_table(labels, label_path)

    high = adata.obs["annotation_confidence"].astype(str).str.lower().eq("high")
    labeled = high & adata.obs["annotation_status"].astype(str).str.lower().eq("labeled")
    labeled_indices = pd.DataFrame(
        {
            "cell_id": adata.obs_names[labeled].astype(str),
            "position": np.flatnonzero(labeled.to_numpy()),
            "label": adata.obs.loc[labeled, "scvi_curated_cell_type"].astype(str).to_numpy(),
        }
    )
    labeled_path = write_table(labeled_indices, table_dir / "scpoli_labeled_indices.tsv")

    cards = []
    label_by_cluster = adata.obs.groupby(cluster_key, observed=True).first()
    for cluster in bundle.clusters:
        row = label_by_cluster.loc[label_by_cluster.index.astype(str) == cluster.cluster_id].iloc[0]
        cards.append(
            {
                "cluster_id": cluster.cluster_id,
                "proposed_label": str(row["scvi_curated_cell_type"]),
                "confidence": str(row["annotation_confidence"]),
                "status": str(row["annotation_status"]),
                "n_cells": cluster.n_cells,
                "n_donors": cluster.n_donors,
                "n_studies": cluster.n_studies,
                "nearest_clusters": ", ".join(cluster.nearest_clusters),
                "top_markers": ", ".join(cluster.top_markers[:20]),
                "dominant_programmes": ", ".join(cluster.dominant_programs),
                "qc_warnings": ", ".join(cluster.qc_warnings),
            }
        )
    cards_frame = pd.DataFrame(cards)
    cards_path = write_table(cards_frame, table_dir / "scvi_annotation_cards.tsv")
    confidence = (
        adata.obs["annotation_confidence"]
        .astype(str)
        .value_counts()
        .rename_axis("confidence")
        .reset_index(name="n_cells")
    )
    overview = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_clusters": int(adata.obs[cluster_key].nunique()),
        "n_labels": int(adata.obs["scvi_curated_cell_type"].nunique()),
        "n_labeled_for_scpoli": int(labeled.sum()),
        "confidence_counts": dict(zip(confidence["confidence"].astype(str), confidence["n_cells"].astype(int))),
        "scvi_embedding_present": "X_scVI" in adata.obsm,
        "cnmf_usage_present": "X_nmf_usage" in adata.obsm,
    }
    sidecar = {
        "schema_version": ATLAS_REPORT_SCHEMA,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "input_h5ad": str(source.resolve()),
        "overview": overview,
        "tables": {
            "subgroup_distribution": str(distribution_path.resolve()),
            "graph_neighbours": str(neighbour_path.resolve()),
            "cnmf_usage": str(programme_path.resolve()),
            "cnmf_program_classification": str(programme_class_path.resolve()),
            "reference_mapping_concordance": str(reference_path.resolve()),
            "annotation_cards": str(cards_path.resolve()),
            "curated_labels": str(label_path.resolve()),
            "scpoli_labeled_indices": str(labeled_path.resolve()),
        },
        "gate": {"status": "blocked", "reason": "user approval is required before scPoli"},
    }
    atomic_json(target.with_suffix(".json"), sidecar)
    sections = [
        "<!doctype html><html><head><meta charset='utf-8'><title>SCOOP scVI annotation report</title>",
        (
            "<style>body{font-family:system-ui,sans-serif;margin:2rem auto;max-width:1200px;line-height:1.45}"
            "table{border-collapse:collapse;width:100%}th,td{border:1px solid #ddd;padding:.4rem;text-align:left}"
            ".gate{border-left:5px solid #2563eb;background:#eff6ff;padding:1rem}</style></head><body>"
        ),
        "<h1>SCOOP scVI annotation report</h1>",
        f"<p>Run: {escape(str(sidecar['created_at']))}; input: {escape(str(source.resolve()))}</p>",
        "<h2>Overview</h2><pre>" + escape(json.dumps(overview, ensure_ascii=False, indent=2)) + "</pre>",
        "<h2>Annotation cards</h2>" + cards_frame.to_html(index=False, escape=True),
        "<h2>Subgroup distribution</h2>" + distributions.head(1000).to_html(index=False, escape=True),
        "<h2>Graph neighbours</h2>" + neighbours.to_html(index=False, escape=True),
        "<h2>cNMF / FastCNMF usage</h2>" + programmes.to_html(index=False, escape=True),
        "<h2>cNMF programme classification</h2>" + programme_classes.to_html(index=False, escape=True),
        "<h2>Reference mapping concordance</h2>" + reference_concordance.to_html(index=False, escape=True),
        "<h2>Confidence and supervision</h2>" + confidence.to_html(index=False, escape=True),
        (
            "<div class='gate'><strong>scPoli blocked.</strong> Review this report and sidecar tables, "
            "then create a hash-bound approval artifact.</div>"
        ),
        "</body></html>",
    ]
    atomic_text(target, "\n".join(sections))
    return target, label_path, labeled_path


def train_scpoli(*, gate_path: str | Path, output_path: str | Path | None = None, n_epochs: int = 100) -> Path:
    verification = verify_approval(gate_path)
    if not verification["ok"]:
        raise PermissionError("scPoli approval gate rejected: " + "; ".join(verification["errors"]))
    payload = verification["payload"]
    artifacts = payload["artifacts"]
    source = Path(artifacts["input_h5ad"]["path"])
    labeled_path = Path(artifacts["labeled_indices"]["path"])
    label_path = Path(artifacts["label_table"]["path"])
    adata = ad.read_h5ad(source)
    labeled_table = _read_labels(labeled_path)
    labels = _read_labels(label_path)
    required_label_columns = {
        "cluster_id",
        "scvi_curated_cell_type",
        "annotation_confidence",
        "annotation_status",
        "annotation_source",
        "annotation_version",
    }
    missing_label_columns = sorted(required_label_columns - set(labels.columns))
    if missing_label_columns:
        raise ValueError("approved label table is incomplete: " + ", ".join(missing_label_columns))
    required_index_columns = {"cell_id", "position", "label"}
    missing_index_columns = sorted(required_index_columns - set(labeled_table.columns))
    if missing_index_columns:
        raise ValueError("labeled indices table is incomplete: " + ", ".join(missing_index_columns))
    labeled_indices = labeled_table["position"].astype(int).to_numpy()
    if len(labeled_indices) == 0:
        raise ValueError("approval contains no high-confidence labeled cells")
    if (labeled_indices < 0).any() or (labeled_indices >= adata.n_obs).any():
        raise IndexError("labeled indices are outside the approved H5AD")
    expected_cell_ids = adata.obs_names[labeled_indices].astype(str).to_numpy()
    if not np.array_equal(expected_cell_ids, labeled_table["cell_id"].astype(str).to_numpy()):
        raise PermissionError("labeled cell IDs do not match their approved H5AD positions")
    confidence = adata.obs.iloc[labeled_indices]["annotation_confidence"].astype(str).str.lower()
    status = adata.obs.iloc[labeled_indices]["annotation_status"].astype(str).str.lower()
    if not confidence.eq("high").all() or not status.eq("labeled").all():
        raise PermissionError("scPoli supervision is restricted to high-confidence labeled cells")
    label_key = "scvi_curated_cell_type"
    approved_labels = adata.obs.iloc[labeled_indices][label_key].astype(str).to_numpy()
    if not np.array_equal(approved_labels, labeled_table["label"].astype(str).to_numpy()):
        raise PermissionError("supervision labels do not match the approved H5AD")
    source_values = adata.obs.iloc[labeled_indices]["annotation_source"].astype(str)
    if not source_values.eq("scvi_cluster_curated").all():
        raise PermissionError("legacy and non-scVI-curated labels cannot supervise scPoli")
    approved_version = str(payload["annotation_version"])
    if not adata.obs.iloc[labeled_indices]["annotation_version"].astype(str).eq(approved_version).all():
        raise PermissionError("annotation version does not match the approval gate")
    for usable_key in ("final_use", "qc_pass"):
        if usable_key in adata.obs and not adata.obs.iloc[labeled_indices][usable_key].astype(bool).all():
            raise PermissionError(f"supervision contains cells rejected by {usable_key}")
    if "doublet_call" in adata.obs:
        doublet = adata.obs.iloc[labeled_indices]["doublet_call"].astype(str).str.lower()
        if doublet.isin({"true", "1", "doublet", "positive"}).any():
            raise PermissionError("supervision contains doublet-called cells")

    cluster_key = next((key for key in ("leiden_scvi", "cluster_identity") if key in adata.obs), None)
    if cluster_key is None:
        raise KeyError("approved H5AD has no supported scVI cluster key")
    if labels["cluster_id"].astype(str).duplicated().any():
        raise ValueError("approved label table has duplicate cluster_id values")
    label_index = labels.assign(cluster_id=labels["cluster_id"].astype(str)).set_index("cluster_id")
    for cluster_id, frame in adata.obs.groupby(cluster_key, observed=True):
        key = str(cluster_id)
        if key not in label_index.index:
            raise PermissionError(f"approved label table is missing cluster {key}")
        record = label_index.loc[key]
        for field in (
            "scvi_curated_cell_type",
            "annotation_confidence",
            "annotation_status",
            "annotation_source",
            "annotation_version",
        ):
            if set(frame[field].astype(str)) != {str(record[field])}:
                raise PermissionError(f"approved label table disagrees with H5AD for cluster {key}: {field}")
    try:
        from scarches.models.scpoli import scPoli
    except ImportError as exc:
        raise RuntimeError("scPoli training requires scarches in the managed SCOOP runtime") from exc

    condition_keys = [key for key in ("sample_id",) if key in adata.obs]
    if not condition_keys:
        raise KeyError("scPoli requires sample_id in obs")
    supervised_labels = set(adata.obs.iloc[labeled_indices][label_key].astype(str))
    unknown = sorted(set(adata.obs[label_key].astype(str)) - supervised_labels)
    original_x = adata.X
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"]
    model = scPoli(
        adata=adata,
        condition_keys=condition_keys,
        cell_type_keys=[label_key],
        labeled_indices=labeled_indices.tolist(),
        unknown_ct_names=unknown,
    )
    model.train(n_epochs=int(n_epochs))
    adata.obsm["X_scPoli"] = model.get_latent(adata, mean=True)
    predictions = model.classify(adata)
    result = predictions[label_key]
    adata.obs["scpoli_predicted_cell_type"] = result["preds"].astype(str)
    adata.obs["scpoli_uncertainty"] = np.asarray(result["uncert"], dtype=float)
    adata.obs["scpoli_weighted_distance"] = np.asarray(result["weighted_distances"], dtype=float)
    adata.X = original_x
    run_root = source.parent.parent
    model_path = run_root / "atlas" / "models" / "scpoli"
    if model_path.exists() and any(model_path.iterdir()):
        raise FileExistsError(f"refusing to overwrite existing scPoli model: {model_path}")
    model_path.parent.mkdir(parents=True, exist_ok=True)
    model.save(str(model_path), overwrite=False)
    target = Path(output_path) if output_path else run_root / "artifacts" / "adata_scpoli.h5ad"
    target.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(target)
    atomic_json(
        run_root / "atlas" / "scpoli-run.json",
        {
            "schema_version": "scoop.scpoli-run.v1",
            "approval_gate": str(Path(gate_path).resolve()),
            "annotation_version": payload["annotation_version"],
            "n_labeled": int(len(labeled_indices)),
            "n_epochs": int(n_epochs),
            "output": str(target.resolve()),
            "model": str(model_path.resolve()),
        },
    )
    return target


def validate_atlas(*, input_path: str | Path, output_path: str | Path) -> Path:
    source = Path(input_path).resolve()
    adata = ad.read_h5ad(source)
    required = ["X_scPoli", "X_scVI"]
    missing = [key for key in required if key not in adata.obsm]
    metrics = {
        "schema_version": "scoop.atlas-validation.v1",
        "input": str(source),
        "ok": not missing,
        "missing_embeddings": missing,
        "n_cells": int(adata.n_obs),
        "n_labels": int(adata.obs["scvi_curated_cell_type"].nunique()) if "scvi_curated_cell_type" in adata.obs else 0,
        "prediction_concordance_on_labeled": None,
    }
    if {"annotation_confidence", "annotation_status", "scvi_curated_cell_type", "scpoli_predicted_cell_type"}.issubset(
        adata.obs
    ):
        mask = adata.obs["annotation_confidence"].astype(str).str.lower().eq("high") & adata.obs[
            "annotation_status"
        ].astype(str).str.lower().eq("labeled")
        if mask.any():
            metrics["prediction_concordance_on_labeled"] = float(
                (
                    adata.obs.loc[mask, "scvi_curated_cell_type"].astype(str)
                    == adata.obs.loc[mask, "scpoli_predicted_cell_type"].astype(str)
                ).mean()
            )
    return atomic_json(Path(output_path), metrics)
