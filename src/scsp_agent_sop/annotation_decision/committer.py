from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Mapping

import pandas as pd
import numpy as np

from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.storage import ensure_dir, register_file, write_table

from .decision_schema import ALLOWED_CONFIDENCE, AnnotationDecision
from .evidence_bundle import AnnotationEvidenceBundle, partition_sha256
from .validator import ValidationResult, validate_annotation_decision


OBS_FIELDS = ("cell_type_lvl1", "cell_type_lvl2", "cell_type_lvl3", "cell_state", "annotation_confidence", "annotation_status",
              "annotation_identity_scope", "broad_annotation_status", "fine_annotation_status")
MERGE_LOG_COLUMNS = ["cluster_id", "cluster_key", "action", "n_cells", "final_label", "parent_unit", "unit_rule"]
CELL_LAYER_FIELDS = {"cell_confidence": "annotation_confidence", "cell_state": "cell_state",
                     "identity_source": "annotation_identity_source"}


def commit_annotation_decisions(
    adata,
    decisions: Iterable[AnnotationDecision],
    evidence_bundle: AnnotationEvidenceBundle,
    *,
    run_root: str | Path,
    output_dir: str | Path,
    min_cells_high_confidence: int = 20,
    min_samples_high_confidence: int = 2,
    tumor_gate: Mapping[str, Any] | None = None,
    cell_layer: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Commit unit decisions; an optional per-cell layer can only lower confidence.

    cell_layer has obs_name plus any of cell_confidence, cell_state and
    identity_source. It applies only to cells whose unit labels are written.
    """
    run_root = Path(run_root)
    output_dir = ensure_dir(output_dir)
    decisions = list(decisions)
    validations: list[ValidationResult] = [
        validate_annotation_decision(decision, evidence_bundle, min_cells_high_confidence=min_cells_high_confidence,
                                     min_samples_high_confidence=min_samples_high_confidence, tumor_gate=tumor_gate)
        for decision in decisions
    ]
    if not adata.obs_names.is_unique:
        raise ValueError("Annotation requires unique observation IDs; preserve library/sample namespaces.")
    if evidence_bundle.partition_sha256 and partition_sha256(adata, evidence_bundle.cluster_key) != evidence_bundle.partition_sha256:
        raise ValueError("Evidence bundle was built on a different partition or obs axis; rerun annotate review.")
    bundle_clusters = evidence_bundle.cluster_by_id()
    units = {}
    positions = []
    covered = np.zeros(adata.n_obs, dtype=bool)
    for decision in decisions:
        if decision.cluster_key not in units:
            if decision.cluster_key not in adata.obs:
                raise KeyError(f"cluster_key {decision.cluster_key!r} is not present in adata.obs")
            labels = adata.obs[decision.cluster_key].astype(str)
            units[decision.cluster_key] = labels.groupby(labels, sort=False, observed=True).indices
        indices = units[decision.cluster_key].get(str(decision.cluster_id))
        if indices is None or len(indices) == 0:
            raise ValueError(f"Decision unit {decision.cluster_key}:{decision.cluster_id} has no observations")
        if covered[indices].any():
            raise ValueError("Annotation decisions overlap; adjudicate one final decision per observation.")
        cluster = bundle_clusters.get(str(decision.cluster_id))
        if decision.cluster_key == evidence_bundle.cluster_key and cluster is not None and cluster.n_cells != len(indices):
            raise ValueError(f"Decision unit {decision.cluster_key}:{decision.cluster_id} has {len(indices)} cells "
                             f"but the evidence bundle has {cluster.n_cells}; rerun annotate review.")
        covered[indices] = True
        positions.append(indices)
    layer = _cell_layer_positions(adata, cell_layer)
    optional = {"malignancy_status": any(d.malignancy_status is not None for d in decisions),
                "malignancy_basis": any(d.malignancy_basis for d in decisions),
                "annotation_technical_kind": any(d.technical_kind is not None for d in decisions),
                "annotation_identity_source": layer is not None and "identity_source" in layer}
    for field in OBS_FIELDS:
        if field not in adata.obs:
            adata.obs[field] = "unassigned" if field != "cell_state" else ""
        elif isinstance(adata.obs[field].dtype, pd.CategoricalDtype):
            adata.obs[field] = adata.obs[field].astype(object)
    for field, declared in optional.items():
        if field in adata.obs:
            adata.obs[field] = adata.obs[field].astype(object)
        elif declared:
            adata.obs[field] = "not_assessed" if field == "malignancy_status" else ""
    optional = [field for field in optional if field in adata.obs]
    if "annotation_reviewed" not in adata.obs:
        adata.obs["annotation_reviewed"] = False
    flags = adata.obs.annotation_reviewed
    if not flags.dropna().isin([True, False, 0, 1]).all():
        raise ValueError("annotation_reviewed must be boolean; missing means not reviewed")
    adata.obs["annotation_reviewed"] = flags.fillna(False).astype(bool)
    columns = {field: adata.obs.columns.get_loc(field) for field in (*OBS_FIELDS, *optional, "annotation_reviewed")}
    decision_rows = []
    validation_rows = []
    merge_rows = []
    unit_confidence = np.full(adata.n_obs, None, dtype=object)
    for decision, validation, indices in zip(decisions, validations, positions, strict=True):
        validation_rows.append(validation.to_dict())
        decision_dict = decision.to_dict()
        decision_dict["accepted"] = validation.accepted
        decision_dict["validation_reasons"] = "; ".join(validation.reasons)
        decision_rows.append(decision_dict)
        technical = validation.accepted and decision.disposition == "technical_exclusion"
        status = "excluded_technical" if technical else "accepted" if validation.accepted else "review_required"
        # A rejected replacement must not leave an old accepted label eligible
        # downstream. Retain its text for provenance, but mark it unresolved.
        adata.obs.iloc[indices, columns["annotation_status"]] = status
        adata.obs.iloc[indices, columns["annotation_reviewed"]] = True
        adata.obs.iloc[indices, columns["annotation_identity_scope"]] = "not_applicable" if technical else decision.identity_scope
        broad_accepted = validation.accepted and decision.identity_scope in {"broad", "fine", "fine_unresolved"}
        fine_accepted = validation.accepted and decision.identity_scope == "fine"
        if technical:
            broad, fine = "not_applicable", "not_applicable"
        else:
            broad = "accepted" if broad_accepted else "review_required"
            fine = ("accepted" if fine_accepted else
                    "unresolved" if validation.accepted and decision.identity_scope == "fine_unresolved" else "review_required")
        adata.obs.iloc[indices, columns["broad_annotation_status"]] = broad
        adata.obs.iloc[indices, columns["fine_annotation_status"]] = fine
        written = validation.accepted or (validation.valid and decision.identity_scope == "unresolved")
        if written:
            values = {"cell_type_lvl1": decision.parent_label, "cell_type_lvl2": decision.canonical_label,
                      "cell_type_lvl3": decision.final_label, "cell_state": decision.cell_state or "",
                      "annotation_confidence": decision.confidence, "malignancy_status": decision.malignancy_status or "not_assessed",
                      "malignancy_basis": ";".join(decision.malignancy_basis), "annotation_technical_kind": decision.technical_kind or "",
                      "annotation_identity_source": "unit"}
            for field, value in values.items():
                if field in columns:
                    adata.obs.iloc[indices, columns[field]] = value
            unit_confidence[indices] = decision.confidence
        action = ("technical_exclusion" if technical else "commit" if validation.accepted
                  else "record_unresolved" if written else "review_required")
        merge_rows.append({"cluster_id": decision.cluster_id, "cluster_key": decision.cluster_key, "action": action,
                           "n_cells": len(indices), "final_label": decision.final_label,
                           "parent_unit": decision.parent_unit, "unit_rule": decision.unit_rule})
        log_decision(
            run_root,
            module="annotation",
            decision=("annotation_technical_exclusion" if technical else
                      "annotation_decision_accepted" if validation.accepted else "annotation_decision_review_required"),
            reason="; ".join(validation.reasons),
            parameters={"cluster_id": decision.cluster_id, "cluster_key": decision.cluster_key, "final_label": decision.final_label},
            evidence={"evidence_refs": decision.evidence_refs},
            human_review_required=validation.review_required,
            review_reason="; ".join(validation.reasons) if validation.review_required else None,
        )
    if layer is not None:
        _apply_cell_layer(adata, layer, unit_confidence, columns)
    decisions_df = pd.DataFrame(decision_rows)
    validation_df = pd.DataFrame(validation_rows)
    merge_df = pd.DataFrame(merge_rows, columns=MERGE_LOG_COLUMNS)
    decisions_path = write_table(decisions_df, output_dir / "annotation_decisions.tsv")
    validation_path = write_table(validation_df, output_dir / "annotation_validation.tsv")
    merge_path = write_table(merge_df, output_dir / "merge_split_log.tsv")
    register_file(adata, key="annotation_decisions", path=decisions_path, schema="annotation_decisions.v1")
    register_file(adata, key="annotation_validation", path=validation_path, schema="annotation_validation.v1")
    register_file(adata, key="annotation_merge_split_log", path=merge_path, schema="annotation_merge_split_log.v1")
    reviewed = adata.obs.annotation_reviewed.fillna(False).astype(bool)
    accepted = adata.obs.annotation_status.isin(["accepted", "labeled"]) & reviewed
    technical = adata.obs.annotation_status.eq("excluded_technical") & reviewed
    adata.uns["annotation_review"] = {
        "scope": "declared_decision_units_not_all_fine_identities",
        "observations": int(adata.n_obs), "reviewed": int(reviewed.sum()),
        "accepted": int(accepted.sum()), "reviewed_unresolved": int((reviewed & ~accepted & ~technical).sum()),
        "technical_excluded": int(technical.sum()),
        "unreviewed": int((~reviewed).sum()), "review_coverage_complete": bool(reviewed.all()),
        "broad_accepted": int((reviewed & adata.obs.broad_annotation_status.eq("accepted")).sum()),
        "fine_accepted": int((reviewed & adata.obs.fine_annotation_status.eq("accepted")).sum()),
    }
    return decisions_df, validation_df, merge_df


def _cell_layer_positions(adata, cell_layer: pd.DataFrame | None) -> pd.DataFrame | None:
    if cell_layer is None:
        return None
    if "obs_name" not in cell_layer or not set(CELL_LAYER_FIELDS) & set(cell_layer.columns):
        raise ValueError(f"cell layer needs obs_name and at least one of {', '.join(CELL_LAYER_FIELDS)}")
    names = cell_layer["obs_name"].astype(str)
    if not names.is_unique:
        raise ValueError("cell layer has duplicate obs_name rows")
    positions = adata.obs_names.get_indexer(names)
    if (positions < 0).any():
        raise ValueError(f"cell layer has {int((positions < 0).sum())} obs_name values absent from the input")
    layer = cell_layer[[field for field in CELL_LAYER_FIELDS if field in cell_layer]].copy()
    layer.index = positions
    if "cell_confidence" in layer:
        values = layer.cell_confidence.dropna().astype(str)
        if not values.isin(ALLOWED_CONFIDENCE).all():
            raise ValueError(f"cell_confidence must be one of {ALLOWED_CONFIDENCE}")
    return layer


def _apply_cell_layer(adata, layer: pd.DataFrame, unit_confidence: np.ndarray, columns: dict) -> None:
    """Apply per-cell values only where this commit wrote the unit's labels."""
    rows = layer[pd.notna(unit_confidence[layer.index.to_numpy()])]
    rank = {value: position for position, value in enumerate(reversed(ALLOWED_CONFIDENCE))}
    for field, column in CELL_LAYER_FIELDS.items():
        if field not in rows:
            continue
        values = rows[field].dropna().astype(str)
        values = values[values.str.strip() != ""]
        if field == "cell_confidence":
            # The unit confidence is a ceiling; a per-cell value can only lower it.
            ceiling = pd.Series(unit_confidence[values.index.to_numpy()], index=values.index)
            values = values.where(values.map(rank) < ceiling.map(rank), ceiling)
        if len(values):
            adata.obs.iloc[values.index.to_numpy(), columns[column]] = values.to_numpy()
