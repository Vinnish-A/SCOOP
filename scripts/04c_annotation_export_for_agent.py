#!/usr/bin/env python
"""Build annotation evidence and a structured decision template.

This script does not call an LLM and does not assign final labels. It prepares
the evidence bundle and blank template that a subagent or analyst fills in.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import anndata as ad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from scsp_agent_sop.annotation_decision.decision_schema import AnnotationDecision
from scsp_agent_sop.annotation_decision.evidence_bundle import build_evidence_bundle
from scsp_agent_sop.config import deep_get, read_yaml, resolve_run_root
from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.storage import ensure_dir, write_json


def _default_input(run_root: Path) -> Path:
    tumor_path = run_root / "artifacts" / "adata_tumor_fastcnv.h5ad"
    if tumor_path.exists():
        return tumor_path
    return run_root / "artifacts" / "adata_annotation_evidence.h5ad"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--cluster-key", default=None)
    ap.add_argument("--output-dir", default=None)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    input_path = Path(args.input) if args.input else _default_input(run_root)
    output_dir = Path(args.output_dir) if args.output_dir else run_root / deep_get(cfg, "annotation.decision.output_dir", "04_annotation/decisions")
    cluster_key = args.cluster_key or deep_get(cfg, "annotation.decision.cluster_key", "cluster_identity")

    adata = ad.read_h5ad(input_path)
    disease_context = deep_get(cfg, "run.disease_context", None)
    is_tumor = bool(deep_get(cfg, "annotation.decision.is_tumor", deep_get(cfg, "tumor_fastcnv.enabled", False)))
    bundle = build_evidence_bundle(
        adata,
        run_id=deep_get(cfg, "run.run_id", run_root.name),
        organism=deep_get(cfg, "run.organism", "human"),
        tissue=deep_get(cfg, "run.tissue", "unknown"),
        disease_context=disease_context,
        is_tumor=is_tumor,
        cluster_key=cluster_key,
    )
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
            ).to_dict()
        )
    ensure_dir(output_dir)
    bundle_path = write_json(bundle.to_dict(), output_dir / "annotation_evidence_bundle.json")
    template_path = write_json({"schema_version": "scoop.annotation_decision_template.v1", "decisions": templates}, output_dir / "annotation_decision_template.json")
    log_decision(
        run_root,
        module="annotation",
        decision="annotation_agent_payload_created",
        reason="Evidence bundle and structured decision template created without assigning labels.",
        parameters={"cluster_key": cluster_key, "input_h5ad": str(input_path)},
        evidence={"evidence_bundle": str(bundle_path), "decision_template": str(template_path), "n_clusters": len(bundle.clusters)},
        human_review_required=True,
        review_reason="Decision template requires subagent- or analyst-edited structured decisions before commit.",
    )


if __name__ == "__main__":
    main()
