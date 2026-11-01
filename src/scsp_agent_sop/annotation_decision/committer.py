from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd
import numpy as np

from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.storage import ensure_dir, register_file, write_table

from .decision_schema import AnnotationDecision
from .evidence_bundle import AnnotationEvidenceBundle
from .validator import ValidationResult, validate_annotation_decision


OBS_FIELDS = ("cell_type_lvl1", "cell_type_lvl2", "cell_type_lvl3", "cell_state", "annotation_confidence", "annotation_status",
              "annotation_identity_scope", "broad_annotation_status", "fine_annotation_status")


def commit_annotation_decisions(
    adata,
    decisions: Iterable[AnnotationDecision],
    evidence_bundle: AnnotationEvidenceBundle,
    *,
    run_root: str | Path,
    output_dir: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    run_root = Path(run_root)
    output_dir = ensure_dir(output_dir)
    decisions = list(decisions)
    validations: list[ValidationResult] = [validate_annotation_decision(decision, evidence_bundle) for decision in decisions]
    if not adata.obs_names.is_unique:
        raise ValueError("Annotation requires unique observation IDs; preserve library/sample namespaces.")
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
        covered[indices] = True
        positions.append(indices)
    for field in OBS_FIELDS:
        if field not in adata.obs:
            adata.obs[field] = "unassigned" if field != "cell_state" else ""
        elif isinstance(adata.obs[field].dtype, pd.CategoricalDtype):
            adata.obs[field] = adata.obs[field].astype(object)
    if "annotation_reviewed" not in adata.obs:
        adata.obs["annotation_reviewed"] = False
    flags = adata.obs.annotation_reviewed
    if not flags.dropna().isin([True, False, 0, 1]).all():
        raise ValueError("annotation_reviewed must be boolean; missing means not reviewed")
    adata.obs["annotation_reviewed"] = flags.fillna(False).astype(bool)
    columns = {field: adata.obs.columns.get_loc(field) for field in (*OBS_FIELDS, "annotation_reviewed")}
    decision_rows = []
    validation_rows = []
    merge_rows = []
    for decision, validation, indices in zip(decisions, validations, positions, strict=True):
        validation_rows.append(validation.to_dict())
        decision_dict = decision.to_dict()
        decision_dict["accepted"] = validation.accepted
        decision_dict["validation_reasons"] = "; ".join(validation.reasons)
        decision_rows.append(decision_dict)
        status = "accepted" if validation.accepted else "review_required"
        # A rejected replacement must not leave an old accepted label eligible
        # downstream. Retain its text for provenance, but mark it unresolved.
        adata.obs.iloc[indices, columns["annotation_status"]] = status
        adata.obs.iloc[indices, columns["annotation_reviewed"]] = True
        adata.obs.iloc[indices, columns["annotation_identity_scope"]] = decision.identity_scope
        broad_accepted = validation.accepted and decision.identity_scope in {"broad", "fine"}
        fine_accepted = validation.accepted and decision.identity_scope == "fine"
        adata.obs.iloc[indices, columns["broad_annotation_status"]] = "accepted" if broad_accepted else "review_required"
        adata.obs.iloc[indices, columns["fine_annotation_status"]] = "accepted" if fine_accepted else "review_required"
        if validation.accepted or (validation.valid and decision.identity_scope == "unresolved"):
            for field, value in {"cell_type_lvl1": decision.parent_label, "cell_type_lvl2": decision.canonical_label,
                                 "cell_type_lvl3": decision.final_label, "cell_state": decision.cell_state or "",
                                 "annotation_confidence": decision.confidence}.items():
                adata.obs.iloc[indices, columns[field]] = value
            merge_rows.append({"cluster_id": decision.cluster_id, "cluster_key": decision.cluster_key,
                               "action": "commit" if validation.accepted else "record_unresolved",
                               "n_cells": len(indices), "final_label": decision.final_label})
        log_decision(
            run_root,
            module="annotation",
            decision="annotation_decision_accepted" if validation.accepted else "annotation_decision_review_required",
            reason="; ".join(validation.reasons),
            parameters={"cluster_id": decision.cluster_id, "cluster_key": decision.cluster_key, "final_label": decision.final_label},
            evidence={"evidence_refs": decision.evidence_refs},
            human_review_required=validation.review_required,
            review_reason="; ".join(validation.reasons) if validation.review_required else None,
        )
    decisions_df = pd.DataFrame(decision_rows)
    validation_df = pd.DataFrame(validation_rows)
    merge_df = pd.DataFrame(merge_rows, columns=["cluster_id", "cluster_key", "action", "n_cells", "final_label"])
    decisions_path = write_table(decisions_df, output_dir / "annotation_decisions.tsv")
    validation_path = write_table(validation_df, output_dir / "annotation_validation.tsv")
    merge_path = write_table(merge_df, output_dir / "merge_split_log.tsv")
    register_file(adata, key="annotation_decisions", path=decisions_path, schema="annotation_decisions.v1")
    register_file(adata, key="annotation_validation", path=validation_path, schema="annotation_validation.v1")
    register_file(adata, key="annotation_merge_split_log", path=merge_path, schema="annotation_merge_split_log.v1")
    reviewed = adata.obs.annotation_reviewed.fillna(False).astype(bool)
    accepted = adata.obs.annotation_status.isin(["accepted", "labeled"]) & reviewed
    adata.uns["annotation_review"] = {
        "scope": "declared_decision_units_not_all_fine_identities",
        "observations": int(adata.n_obs), "reviewed": int(reviewed.sum()),
        "accepted": int(accepted.sum()), "reviewed_unresolved": int((reviewed & ~accepted).sum()),
        "unreviewed": int((~reviewed).sum()), "review_coverage_complete": bool(reviewed.all()),
        "broad_accepted": int((reviewed & adata.obs.broad_annotation_status.eq("accepted")).sum()),
        "fine_accepted": int((reviewed & adata.obs.fine_annotation_status.eq("accepted")).sum()),
    }
    return decisions_df, validation_df, merge_df
