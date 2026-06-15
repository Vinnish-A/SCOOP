from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd

from scsp_agent_sop.decision_log import log_decision
from scsp_agent_sop.storage import ensure_dir, register_file, write_table

from .decision_schema import AnnotationDecision
from .evidence_bundle import AnnotationEvidenceBundle
from .validator import ValidationResult, validate_annotation_decision


OBS_FIELDS = ("cell_type_lvl1", "cell_type_lvl2", "cell_type_lvl3", "cell_state", "annotation_confidence", "annotation_status")


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
    for field in OBS_FIELDS:
        if field not in adata.obs:
            adata.obs[field] = "unassigned" if field != "cell_state" else ""
    decision_rows = []
    validation_rows = []
    merge_rows = []
    for decision, validation in zip(decisions, validations, strict=True):
        validation_rows.append(validation.to_dict())
        decision_dict = decision.to_dict()
        decision_dict["accepted"] = validation.accepted
        decision_dict["validation_reasons"] = "; ".join(validation.reasons)
        decision_rows.append(decision_dict)
        status = "accepted" if validation.accepted else "review_required"
        if validation.accepted:
            if decision.cluster_key not in adata.obs:
                raise KeyError(f"cluster_key {decision.cluster_key!r} is not present in adata.obs")
            mask = adata.obs[decision.cluster_key].astype(str) == str(decision.cluster_id)
            adata.obs.loc[mask, "cell_type_lvl1"] = decision.parent_label
            adata.obs.loc[mask, "cell_type_lvl2"] = decision.canonical_label
            adata.obs.loc[mask, "cell_type_lvl3"] = decision.final_label
            adata.obs.loc[mask, "cell_state"] = decision.cell_state or ""
            adata.obs.loc[mask, "annotation_confidence"] = decision.confidence
            adata.obs.loc[mask, "annotation_status"] = status
            merge_rows.append({"cluster_id": decision.cluster_id, "cluster_key": decision.cluster_key, "action": "commit", "n_cells": int(mask.sum()), "final_label": decision.final_label})
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
    return decisions_df, validation_df, merge_df
