#!/usr/bin/env python
"""Validate and commit structured annotation decisions into H5AD obs."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import anndata as ad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from scsp_agent_sop.annotation_decision.committer import commit_annotation_decisions
from scsp_agent_sop.annotation_decision.decision_schema import decisions_from_json
from scsp_agent_sop.annotation_decision.evidence_bundle import AnnotationEvidenceBundle, build_evidence_bundle
from scsp_agent_sop.config import deep_get, read_yaml, resolve_run_root
from scsp_agent_sop.storage import ensure_dir, init_file_registry


def _default_input(run_root: Path) -> Path:
    tumor_path = run_root / "artifacts" / "adata_tumor_fastcnv.h5ad"
    if tumor_path.exists():
        return tumor_path
    return run_root / "artifacts" / "adata_annotation_evidence.h5ad"


def _load_or_build_bundle(path: Path | None, adata, cfg: dict, run_root: Path, cluster_key: str) -> AnnotationEvidenceBundle:
    if path is not None and path.exists():
        return AnnotationEvidenceBundle.from_json(path.read_text(encoding="utf-8"))
    return build_evidence_bundle(
        adata,
        run_id=deep_get(cfg, "run.run_id", run_root.name),
        organism=deep_get(cfg, "run.organism", "human"),
        tissue=deep_get(cfg, "run.tissue", "unknown"),
        disease_context=deep_get(cfg, "run.disease_context", None),
        is_tumor=bool(deep_get(cfg, "annotation.decision.is_tumor", deep_get(cfg, "tumor_fastcnv.enabled", False))),
        cluster_key=cluster_key,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--decisions", required=True)
    ap.add_argument("--evidence-bundle", default=None)
    ap.add_argument("--output", default=None)
    ap.add_argument("--cluster-key", default=None)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    input_path = Path(args.input) if args.input else _default_input(run_root)
    decisions_path = Path(args.decisions)
    output_path = Path(args.output) if args.output else run_root / "artifacts" / "adata_annotation_committed.h5ad"
    output_dir = run_root / deep_get(cfg, "annotation.decision.output_dir", "04_annotation/decisions")
    cluster_key = args.cluster_key or deep_get(cfg, "annotation.decision.cluster_key", "cluster_identity")

    adata = ad.read_h5ad(input_path)
    init_file_registry(adata, deep_get(cfg, "run.run_id", run_root.name))
    bundle_path = Path(args.evidence_bundle) if args.evidence_bundle else decisions_path.parent / "annotation_evidence_bundle.json"
    evidence_bundle = _load_or_build_bundle(bundle_path, adata, cfg, run_root, cluster_key)
    decisions = decisions_from_json(decisions_path.read_text(encoding="utf-8"))
    commit_annotation_decisions(adata, decisions, evidence_bundle, run_root=run_root, output_dir=output_dir)

    ensure_dir(output_path.parent)
    adata.write_h5ad(output_path)


if __name__ == "__main__":
    main()
