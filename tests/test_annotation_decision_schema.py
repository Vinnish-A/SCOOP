from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scsp_agent_sop.annotation_decision.decision_schema import AnnotationDecision
from scsp_agent_sop.annotation_decision.evidence_bundle import AnnotationEvidenceBundle, build_evidence_bundle


def test_evidence_bundle_builds_with_missing_optional_fields() -> None:
    adata = AnnData(np.ones((4, 3)))
    adata.obs = pd.DataFrame({"cluster_identity": ["0", "0", "1", "1"], "sample_id": ["s1", "s2", "s1", "s2"]}, index=[f"c{i}" for i in range(4)])
    bundle = build_evidence_bundle(
        adata,
        run_id="run1",
        organism="human",
        tissue="brain",
        disease_context=None,
        is_tumor=False,
        cluster_key="cluster_identity",
    )
    assert isinstance(bundle, AnnotationEvidenceBundle)
    assert len(bundle.clusters) == 2
    assert bundle.clusters[0].marker_refs == ()


def test_evidence_bundle_includes_fastcnvpy_fields_when_present() -> None:
    adata = AnnData(np.ones((3, 2)))
    adata.obs = pd.DataFrame(
        {
            "cluster_identity": ["0", "0", "1"],
            "sample_id": ["s1", "s2", "s1"],
            "fastcnv_tumor_evidence": ["tumor", "tumor", "normal"],
            "fastcnv_cnv_fraction": [0.4, 0.5, 0.01],
        },
        index=[f"c{i}" for i in range(3)],
    )
    bundle = build_evidence_bundle(
        adata,
        run_id="run1",
        organism="human",
        tissue="brain",
        disease_context="glioblastoma",
        is_tumor=True,
        cluster_key="cluster_identity",
    )
    cluster0 = bundle.cluster_by_id()["0"]
    assert "fastcnv_tumor_evidence" in cluster0.cnv_refs
    assert cluster0.cnv_summary["has_tumor_evidence"] is True


def test_invalid_confidence_is_rejected() -> None:
    with pytest.raises(ValueError):
        AnnotationDecision.from_dict(
            {
                "schema_version": "scoop.annotation_decision.v1",
                "run_id": "run1",
                "cluster_id": "0",
                "cluster_key": "cluster_identity",
                "parent_label": "Immune",
                "canonical_label": "T cell",
                "cell_state": None,
                "functional_modifier": None,
                "final_label": "T cell",
                "confidence": "certain",
                "evidence_refs": {},
                "positive_markers": [],
                "negative_markers_absent": [],
                "conflicts": [],
                "review_required": False,
                "reason": "bad confidence",
            }
        )
