#!/usr/bin/env python
"""Build annotation evidence and a structured decision template.

This script does not call an LLM and does not assign final labels. It prepares
the evidence bundle and blank template that a subagent or analyst fills in.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd

from scsp_agent_sop.annotation_decision.decision_schema import AnnotationDecision
from scsp_agent_sop.annotation_decision.evidence_bundle import build_evidence_bundle, bound_evidence_files
from scsp_agent_sop.config import ObsKeys, deep_get, read_yaml, resolve_obs_keys, resolve_run_root
from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.storage import ensure_dir, sha256_file, write_json

from scoop.results import archive_receipt


def _default_input(run_root: Path) -> Path:
    tumor_path = run_root / "artifacts" / "adata_tumor_fastcnv.h5ad"
    if tumor_path.exists():
        return tumor_path
    return run_root / "artifacts" / "adata_annotation_evidence.h5ad"


def composition_keys(cfg: dict, obs) -> ObsKeys:
    """annotation.decision.{sample,donor,study}_key override run.obs_keys, then the legacy keys.*."""
    return resolve_obs_keys(cfg, obs, overrides={role: f"annotation.decision.{role}_key" for role in ("sample", "donor", "study")})


def join_cluster_table(adata, path: Path, cluster_key: str) -> None:
    """Take the decision partition from a table keyed by obs_name, not an H5AD column (C1-08)."""
    suffix = path.suffix.lower()
    table = pd.read_parquet(path) if suffix == ".parquet" else pd.read_csv(path, sep="," if suffix == ".csv" else "\t", dtype=str)
    if "obs_name" not in table or cluster_key not in table:
        raise KeyError(f"{path} needs obs_name and {cluster_key!r} columns")
    if cluster_key in adata.obs:
        raise ValueError(f"{cluster_key!r} is both an obs column and a --cluster-table column; keep one source")
    if table.obs_name.duplicated().any():
        raise ValueError(f"{path} repeats obs_name values")
    labels = table.set_index(table.obs_name.astype(str))[cluster_key].reindex(adata.obs_names.astype(str))
    if labels.isna().any():
        raise ValueError(f"{path} gives no {cluster_key!r} for {int(labels.isna().sum())} of {adata.n_obs} observations")
    adata.obs[cluster_key] = pd.Categorical(labels.astype(str).to_numpy())


def bundle_options(cfg: dict, keys: ObsKeys) -> dict:
    """Composition keys and per-sample support threshold for the evidence bundle."""
    return {"sample_key": keys.sample, "donor_key": keys.donor, "study_key": keys.study,
            "min_cells_per_sample": int(deep_get(cfg, "annotation.decision.min_cells_per_sample", 20))}


def main() -> dict:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--cluster-key", default=None)
    ap.add_argument("--cluster-table", default=None, help="Table with obs_name and the cluster key column (C1-08)")
    ap.add_argument("--output-dir", default=None)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    input_path = Path(args.input) if args.input else _default_input(run_root)
    output_dir = (
        Path(args.output_dir)
        if args.output_dir
        else run_root / deep_get(cfg, "annotation.decision.output_dir", "annotation/decisions")
    )
    cluster_key = args.cluster_key or deep_get(cfg, "annotation.decision.cluster_key", "cluster_identity")

    adata = ad.read_h5ad(input_path)
    if args.cluster_table:
        join_cluster_table(adata, Path(args.cluster_table), cluster_key)
    disease_context = deep_get(cfg, "run.disease_context", None)
    is_tumor = bool(deep_get(cfg, "annotation.decision.is_tumor", deep_get(cfg, "tumor_fastcnv.enabled", False)))
    keys = composition_keys(cfg, adata.obs)
    bundle = build_evidence_bundle(
        adata,
        run_id=deep_get(cfg, "run.run_id", run_root.name),
        organism=deep_get(cfg, "run.organism", "human"),
        tissue=deep_get(cfg, "run.tissue", "unknown"),
        disease_context=disease_context,
        is_tumor=is_tumor,
        cluster_key=cluster_key,
        **bundle_options(cfg, keys),
    )
    cited = {ref for cluster in bundle.clusters for refs in (cluster.deg_refs, cluster.marker_refs, cluster.nmf_refs,
                                                             cluster.reference_refs, cluster.cnv_refs) for ref in refs}
    evidence_inputs = bound_evidence_files(bundle, cited)
    ensure_dir(output_dir)
    bundle_file = output_dir / "annotation_evidence_bundle.json"
    template_file = output_dir / "annotation_decision_template.json"
    # Decisions drafted against an earlier bundle stay checkable (C1-09).
    archive_receipt(bundle_file)
    archive_receipt(template_file)
    bundle_path = write_json(bundle.to_dict(), bundle_file)
    bundle_sha256 = sha256_file(bundle_path)
    templates = []
    for cluster in bundle.clusters:
        templates.append(
            AnnotationDecision(
                schema_version="scoop.annotation_decision.v1",
                run_id=bundle.run_id,
                cluster_id=cluster.cluster_id,
                cluster_key=cluster.cluster_key,
                parent_label="review_required",
                canonical_label="review_required",
                cell_state=None,
                functional_modifier=None,
                final_label="review_required",
                confidence="low",
                evidence_refs={
                    "deg": cluster.deg_refs,
                    "markers": cluster.marker_refs,
                    "programs": cluster.nmf_refs,
                    "references": cluster.reference_refs,
                    "cnv": cluster.cnv_refs,
                },
                positive_markers=cluster.top_markers,
                negative_markers_absent=(),
                conflicts=(),
                review_required=True,
                reason="template only; subagent or analyst decision required",
                bundle_sha256=bundle_sha256,
            ).to_dict()
        )
    template_path = write_json(
        {"schema_version": "scoop.annotation_decision_template.v1", "decisions": templates},
        template_file,
    )
    without_deg = [cluster.cluster_id for cluster in bundle.clusters if not cluster.deg_refs]
    log_decision(
        run_root,
        module="annotation",
        decision="annotation_agent_payload_created",
        reason="Evidence bundle and structured decision template created without assigning labels.",
        parameters={"cluster_key": cluster_key, "input_h5ad": str(input_path), "obs_keys": keys.record()},
        evidence={
            "evidence_bundle": str(bundle_path),
            "decision_template": str(template_path),
            "n_clusters": len(bundle.clusters),
        },
        human_review_required=True,
        review_reason="Decision template requires subagent- or analyst-edited structured decisions before commit.",
    )
    return {"inputs": [str(input_path.resolve()), *([str(Path(args.cluster_table).resolve())] if args.cluster_table else []),
                       *map(str, evidence_inputs)],
            "outputs": [str(bundle_path.resolve()), str(template_path.resolve())],
            "execution": {"status": "completed"},
            "validation": {"status": "limited", "scope": "evidence_bundle_and_template_preparation"},
            "interpretation": {"status": "not_assessed"},
            "obs_keys": keys.record(),
            "coverage": [{"scope": "prepared_review_units_not_decisions", "unit": "numeric_clusters",
                          "total": len(bundle.clusters), "processed": len(templates)},
                         {"scope": "units_with_partition_deg", "unit": "numeric_clusters",
                          "total": len(bundle.clusters), "processed": len(bundle.clusters) - len(without_deg)}],
            "limitations": ["Templates are prepared; biological identities have not been adjudicated by this command.",
                            *([f"{len(without_deg)} units have no registered DEG table (label_key {cluster_key}) testing their group; "
                               "cite another partition's group as ref::group or record the gap."] if without_deg else [])]}


if __name__ == "__main__":
    main()
