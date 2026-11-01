#!/usr/bin/env python
"""Validate and commit structured annotation decisions into H5AD obs."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd

from scsp_agent_sop.annotation_decision.committer import OBS_FIELDS, commit_annotation_decisions
from scsp_agent_sop.annotation_decision.decision_schema import decisions_from_json
from scsp_agent_sop.annotation_decision.evidence_bundle import (
    AnnotationEvidenceBundle,
    bind_legacy_bundle,
    bound_evidence_files,
    build_evidence_bundle,
)
from scsp_agent_sop.config import deep_get, project_root_from_run_root, read_yaml, resolve_run_root
from scsp_agent_sop.storage import (changed_obs_columns, ensure_dir, init_file_registry, register_file, sha256_file,
                                    write_h5ad_atomic, write_obs_sidecar)

from scoop.commands.annotation_review import bundle_options, composition_keys, join_cluster_table


def _default_input(run_root: Path) -> Path:
    tumor_path = run_root / "artifacts" / "adata_tumor_fastcnv.h5ad"
    if tumor_path.exists():
        return tumor_path
    return run_root / "artifacts" / "adata_annotation_evidence.h5ad"


def _load_or_build_bundle(
    path: Path | None, adata, cfg: dict, run_root: Path, cluster_key: str
) -> AnnotationEvidenceBundle:
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
        **bundle_options(cfg, composition_keys(cfg, adata.obs)),
    )


def _tumor_gate(cfg: dict, run_root: Path) -> dict:
    schema_path = project_root_from_run_root(run_root) / deep_get(cfg, "annotation.decision.schema", "configs/annotation_decision_schema.yaml")
    if not schema_path.is_file():
        schema_path = Path(__file__).resolve().parents[1] / "data" / "annotation-decision-schema.yaml"
    return dict(read_yaml(schema_path).get("tumor_gate") or {})


def _read_cell_layer(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, sep="," if suffix == ".csv" else "\t", dtype={"obs_name": str})


def main() -> dict:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input", default=None)
    ap.add_argument("--decisions", required=True)
    ap.add_argument("--evidence-bundle", default=None)
    ap.add_argument("--output", default=None)
    ap.add_argument("--cluster-key", default=None)
    ap.add_argument("--cluster-table", default=None,
                    help="Table with obs_name and the cluster key column, for an input without that obs column (C1-08)")
    ap.add_argument("--output-format", choices=["h5ad", "obs-parquet"], default="h5ad",
                    help="obs-parquet writes only the committed and changed obs columns keyed by obs_name, plus a JSON receipt")
    ap.add_argument("--cell-layer", default=None,
                    help="Table with obs_name and cell_confidence/cell_state/identity_source; confidence can only be lowered")
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    run_root = resolve_run_root(args.config, cfg)
    input_path = Path(args.input) if args.input else _default_input(run_root)
    decisions_path = Path(args.decisions)
    sidecar = args.output_format == "obs-parquet"
    default_output = "annotation_committed_obs.parquet" if sidecar else "adata_annotation_committed.h5ad"
    output_path = Path(args.output) if args.output else run_root / "artifacts" / default_output
    output_dir = run_root / deep_get(cfg, "annotation.decision.output_dir", "annotation/decisions")
    cluster_key = args.cluster_key or deep_get(cfg, "annotation.decision.cluster_key", "cluster_identity")

    adata = ad.read_h5ad(input_path)
    init_file_registry(adata, deep_get(cfg, "run.run_id", run_root.name))
    if args.cluster_table:
        join_cluster_table(adata, Path(args.cluster_table), cluster_key)
    input_obs = adata.obs.copy() if sidecar else None
    bundle_path = (
        Path(args.evidence_bundle)
        if args.evidence_bundle
        else decisions_path.parent / "annotation_evidence_bundle.json"
    )
    evidence_bundle = _load_or_build_bundle(bundle_path, adata, cfg, run_root, cluster_key)
    decisions = decisions_from_json(decisions_path.read_text(encoding="utf-8"))
    bundle_sha256 = sha256_file(bundle_path) if bundle_path.is_file() else None
    stale = sorted({d.cluster_id for d in decisions if d.bundle_sha256 and d.bundle_sha256 != bundle_sha256})
    if stale:
        raise ValueError(f"Decisions for units {', '.join(stale[:10])} were drafted against a different evidence bundle "
                         f"than {bundle_path}; re-draft them or restore that bundle from history/.")
    # v1 bundles predate partition and DEG-content binding; rebuild both from
    # this input, so only the per-unit cell counts tie them to the review.
    legacy = evidence_bundle.partition_sha256 is None
    if legacy:
        evidence_bundle = bind_legacy_bundle(evidence_bundle, adata)
    inputs = [str(input_path.resolve()), str(decisions_path.resolve()),
              *([str(Path(args.cluster_table).resolve())] if args.cluster_table else [])]
    if bundle_path.is_file():
        inputs.append(str(bundle_path.resolve()))
    cell_layer = None
    if args.cell_layer:
        cell_layer_path = Path(args.cell_layer)
        cell_layer = _read_cell_layer(cell_layer_path)
        inputs.append(str(cell_layer_path.resolve()))
    cited = {ref for decision in decisions for refs in decision.evidence_refs.values() for ref in refs}
    inputs.extend(str(path) for path in bound_evidence_files(evidence_bundle, cited))
    commit_annotation_decisions(
        adata, decisions, evidence_bundle, run_root=run_root, output_dir=output_dir,
        min_cells_high_confidence=int(deep_get(cfg, "annotation.decision.min_cells_high_confidence", 20)),
        min_samples_high_confidence=int(deep_get(cfg, "annotation.decision.min_samples_high_confidence", 2)),
        tumor_gate=_tumor_gate(cfg, run_root), cell_layer=cell_layer,
    )
    if cell_layer is not None:
        register_file(adata, key="annotation_cell_layer", path=cell_layer_path.resolve(), schema="annotation_cell_layer.v1")
    ensure_dir(output_path.parent)
    outputs = [str(output_path.resolve())]
    if sidecar:
        # the committed label state as a whole, plus optional fields the commit declared
        changed = changed_obs_columns(input_obs, adata)
        columns = [c for c in adata.obs.columns if c in OBS_FIELDS or c == "annotation_reviewed" or c in changed]
        outputs.append(str(write_obs_sidecar(adata, columns, output_path, key="annotation_committed_obs").resolve()))
    else:
        write_h5ad_atomic(adata, output_path)
    coverage = dict(adata.uns["annotation_review"])
    checks = [{"name": "decision_to_observation_mapping", "status": "passed"},
              {"name": "bundle_partition_binding", "status": "limited" if legacy else "passed"},
              {"name": "decision_bundle_hash", "status": "passed" if all(d.bundle_sha256 for d in decisions)
               else "limited" if any(d.bundle_sha256 for d in decisions) else "not_assessed"}]
    return {"inputs": inputs,
            "outputs": [*outputs, str((output_dir / "annotation_decisions.tsv").resolve()),
                        str((output_dir / "annotation_validation.tsv").resolve()),
                        str((output_dir / "merge_split_log.tsv").resolve())],
            "execution": {"status": "completed"},
            "validation": {"status": "passed" if coverage["unreviewed"] == 0 and not legacy else "limited",
                           "scope": "decision_schema_and_observation_review_coverage",
                           "checks": checks},
            "interpretation": {"status": "not_assessed"},
            "coverage": [{"scope": "annotation_review", "unit": "observation_views",
                          "total": int(adata.n_obs), "processed": int(coverage["reviewed"])}],
            "review_counts": coverage,
            "output_format": args.output_format,
            "limitations": ["Review coverage includes explicit unresolved decisions and does not establish all fine identities.",
                            *(["Legacy v1 evidence bundle: partition and DEG coverage were rebuilt from this input, "
                               "not bound at review time."] if legacy else [])]}


if __name__ == "__main__":
    main()
