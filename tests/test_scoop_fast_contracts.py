from __future__ import annotations

from scoop_fast import ArtifactBundle
from scoop_fast.registry import get_engine, list_engines


def test_fast_registry_exposes_expected_engine_ids() -> None:
    ids = {engine.engine_id for engine in list_engines()}
    assert {
        "fastde.markers",
        "fastde.pseudobulk_deseq2",
        "fastde.abundance",
        "fastcnmf.programs",
        "fastcnvpy.tumor_pooled",
        "fastcopykat.cnv_prediction",
    }.issubset(ids)
    assert get_engine("fastde.markers").task_type == "marker_genes"


def test_artifact_bundle_json_round_trip() -> None:
    bundle = ArtifactBundle(
        schema_version="scoop.artifact_bundle.v1",
        engine_id="fastde.markers",
        task_type="marker_genes",
        run_id="run1",
        status="completed",
        inputs={"h5ad": "input.h5ad"},
        outputs={"markers": "markers.tsv"},
        quality={"ok": True},
        timings={"seconds": 1.2},
        registry_patch={"tables": {"markers": {"path": "markers.tsv"}}},
        decision_log_patch={"decision": "markers_completed"},
    )
    restored = ArtifactBundle.from_json(bundle.to_json())
    assert restored == bundle
