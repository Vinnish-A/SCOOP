import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scoop.commands.deliver import review_coverage
from scsp_agent_sop.manifest import build_manifest
import sys
import yaml
from scoop.commands import deliver


def test_explicit_unresolved_review_is_deliverable_but_unreviewed_is_not():
    obj = AnnData(np.zeros((2, 2)), obs=pd.DataFrame({"annotation_reviewed": [True, True],
        "annotation_status": ["accepted", "review_required"], "cell_type_lvl3": ["Mast", "unassigned"]}, index=["a", "b"]))
    assert review_coverage(obj)["reviewed_unresolved"] == 1
    obj.obs.loc["b", "annotation_reviewed"] = False
    with pytest.raises(ValueError, match="Unreviewed"):
        review_coverage(obj)


def test_accepted_status_alone_does_not_prove_label_assignment():
    obj = AnnData(np.zeros((1, 1)), obs=pd.DataFrame({"annotation_reviewed": [True],
        "annotation_status": ["accepted"], "cell_type_lvl3": ["unassigned"]}, index=["a"]))
    with pytest.raises(ValueError, match="assigned labels"):
        review_coverage(obj)


def test_manifest_uses_registered_artifacts_and_reports_missing(tmp_path):
    included, cache = tmp_path / "final.h5ad", tmp_path / "obsolete_cache"
    included.write_bytes(b"delivery")
    cache.write_bytes(b"stale")
    manifest = build_manifest(tmp_path, paths=[included, tmp_path / "missing"])
    assert [x["path"] for x in manifest["files"]] == ["final.h5ad"]
    assert manifest["missing_registered_artifacts"] == [str(tmp_path / "missing")]


def test_default_delivery_never_falls_back_to_existing_qc(tmp_path, monkeypatch):
    run = tmp_path / "runs/test"
    (run / "config").mkdir(parents=True)
    (run / "artifacts").mkdir()
    config = run / "config/run.yaml"
    config.write_text(yaml.safe_dump({"run": {"run_id": "test"}}))
    AnnData(np.zeros((1, 1))).write_h5ad(run / "artifacts/adata_qc.h5ad")
    monkeypatch.setattr(sys, "argv", ["deliver", "--config", str(config)])
    with pytest.raises(FileNotFoundError, match="reviewed-annotation"):
        deliver.main()
    assert not (run / "artifacts/final_adata.h5ad").exists()


def test_reviewed_unresolved_delivery_roundtrips_with_scoped_receipt(tmp_path, monkeypatch):
    run = tmp_path / "runs/test"
    (run / "config").mkdir(parents=True)
    (run / "artifacts").mkdir()
    config = run / "config/run.yaml"
    config.write_text(yaml.safe_dump({"run": {"run_id": "test"}}))
    obj = AnnData(np.zeros((2, 2)), obs=pd.DataFrame({"annotation_reviewed": [True, True],
        "annotation_status": ["accepted", "review_required"], "cell_type_lvl3": ["Mast", "Glial_Unresolved"]}, index=["a", "b"]))
    obj.write_h5ad(run / "artifacts/adata_annotation_committed.h5ad")
    monkeypatch.setattr(sys, "argv", ["deliver", "--config", str(config)])
    result = deliver.main()
    from anndata import read_h5ad
    saved = read_h5ad(run / "artifacts/final_adata.h5ad")
    assert saved.uns["delivery_scope"]["reviewed_unresolved"] == 1
    assert not saved.uns["delivery_scope"]["all_fine_identities_resolved"]
    assert result["validation"]["status"] == "passed"
    assert result["interpretation"]["status"] == "not_assessed"


def test_reviewed_delivery_reads_the_obs_sidecar_commit_declared(tmp_path, monkeypatch):
    import json
    from anndata import read_h5ad
    from scsp_agent_sop.storage import write_obs_sidecar

    run = tmp_path / "runs/test"
    (run / "config").mkdir(parents=True)
    (run / "artifacts").mkdir()
    (run / "state/stages").mkdir(parents=True)
    config = run / "config/run.yaml"
    config.write_text(yaml.safe_dump({"run": {"run_id": "test"}}))
    base = run / "artifacts/adata_annotation_evidence.h5ad"
    obj = AnnData(np.zeros((2, 2)), obs=pd.DataFrame(index=["a", "b"]))
    obj.write_h5ad(base)
    # an older full-object commit that the sidecar commit superseded
    stale = obj.copy()
    stale.obs["annotation_reviewed"], stale.obs["annotation_status"], stale.obs["cell_type_lvl3"] = False, "unassigned", ""
    stale.write_h5ad(run / "artifacts/adata_annotation_committed.h5ad")
    obj.obs["annotation_reviewed"] = [True, True]
    obj.obs["annotation_status"] = ["accepted", "review_required"]
    obj.obs["cell_type_lvl3"] = ["Mast", "Glial_Unresolved"]
    sidecar = run / "artifacts/annotation_committed_obs.parquet"
    write_obs_sidecar(obj, ["annotation_reviewed", "annotation_status", "cell_type_lvl3"], sidecar, base=base,
                      key="annotation_committed_obs")
    (run / "state/stages/annotation-commit.json").write_text(json.dumps(
        {"stage": "annotation-commit", "status": "committed", "outputs": [{"path": str(sidecar)}]}))
    monkeypatch.setattr(sys, "argv", ["deliver", "--config", str(config)])
    deliver.main()
    saved = read_h5ad(run / "artifacts/final_adata.h5ad")
    assert saved.obs.cell_type_lvl3.tolist() == ["Mast", "Glial_Unresolved"]
    assert saved.uns["delivery_scope"]["input"] == str(sidecar.resolve())
