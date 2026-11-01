from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import yaml

from scoop.workflow import run_command_module
from scsp_agent_sop.decision_log import log_decision


def test_qc_receipt_retains_review_and_h5ad_output(tmp_path, monkeypatch):
    run = tmp_path / "runs/test"
    (run / "config").mkdir(parents=True)
    config = run / "config/run.yaml"
    config.write_text(yaml.safe_dump({"run": {"run_id": "test"}}))

    def main():
        (run / "artifacts").mkdir()
        (run / "artifacts/adata_qc.h5ad").write_bytes(b"test")
        log_decision(
            run,
            module="qc",
            decision="qc_complete",
            reason="filtered",
            human_review_required=True,
            review_reason="Sample loses 40% of observations",
        )

    monkeypatch.setitem(sys.modules, "receipt_test", SimpleNamespace(main=main))
    assert run_command_module("qc", "receipt_test", ["--config", str(config)]) == 0
    receipt = json.loads((run / "state/stages/qc.json").read_text())
    assert receipt["status"] == "awaiting_review"
    assert "40%" in receipt["review_reason"]
    assert any(record["path"].endswith("adata_qc.h5ad") for record in receipt["outputs"])


def test_spatial_preparation_is_not_analysis_completion(tmp_path, monkeypatch):
    run = tmp_path / "runs/test"
    (run / "config").mkdir(parents=True)
    config = run / "config/run.yaml"
    config.write_text(yaml.safe_dump({"spatial": {"enabled": True}}))
    monkeypatch.setitem(sys.modules, "receipt_test", SimpleNamespace(main=lambda: None))
    run_command_module("spatial", "receipt_test", ["--config", str(config), "--dry-run"])
    assert json.loads((run / "state/stages/spatial.json").read_text())["status"] == "prepared"


def test_stage_reporting_not_run_is_recorded_as_skipped(tmp_path, monkeypatch):
    run = tmp_path / "runs/test"
    (run / "config").mkdir(parents=True)
    config = run / "config/run.yaml"
    config.write_text(yaml.safe_dump({"run": {"run_id": "test"}}))
    result = {"execution": {"status": "not_run", "reason": "no eligible reference"}}
    monkeypatch.setitem(sys.modules, "receipt_test", SimpleNamespace(main=lambda: result))
    assert run_command_module("tumor_cnv", "receipt_test", ["--config", str(config)]) == 0
    receipt = json.loads((run / "state/stages/tumor_cnv.json").read_text())
    assert receipt["status"] == "skipped"
    assert receipt["execution"]["status"] == "not_run"
    assert receipt["execution"]["reason"] == "no eligible reference"


def test_receipts_go_stale_on_the_config_sections_and_data_their_stage_reads(tmp_path, monkeypatch):
    from scoop.results import inspect_run

    run = tmp_path / "runs/test"
    (run / "config").mkdir(parents=True)
    (run / "artifacts").mkdir()
    config = run / "config/run.yaml"
    config.write_text(yaml.safe_dump({"run": {"run_id": "test"}, "qc": {"min_genes": 200}, "spatial": {"enabled": False}}))
    (run / "artifacts/validated_adata.h5ad").write_bytes(b"validated")

    def main():  # like qc, returns nothing and reads its default input
        (run / "artifacts/adata_qc.h5ad").write_bytes(b"qc")

    monkeypatch.setitem(sys.modules, "receipt_test", SimpleNamespace(main=main))
    assert run_command_module("qc", "receipt_test", ["--config", str(config)]) == 0
    receipt = json.loads((run / "state/stages/qc.json").read_text())
    data = [record for record in receipt["inputs"] if record["role"] == "data"]
    assert [Path(record["path"]).name for record in data] == ["validated_adata.h5ad"]
    configuration = receipt["inputs"][0]
    assert configuration["config_sections"] == ["run", "keys", "qc"] and configuration["sha256"] == receipt["config_hash"]

    edited = yaml.safe_load(config.read_text())
    edited["spatial"]["enabled"] = True
    edited["ccc"] = {"enabled": False}
    config.write_text(yaml.safe_dump(edited))
    assert inspect_run(run, verify=True)["stale_or_missing_stages"] == []
    edited["qc"]["min_genes"] = 100
    config.write_text(yaml.safe_dump(edited))
    assert inspect_run(run, verify=True)["stale_or_missing_stages"] == ["qc"]
    config.write_text(yaml.safe_dump({**edited, "qc": {"min_genes": 200}}))
    (run / "artifacts/validated_adata.h5ad").write_bytes(b"replaced")
    assert inspect_run(run, verify=True)["stale_or_missing_stages"] == ["qc"]
    assert run_command_module("qc", "receipt_test", ["--config", str(config)]) == 0
    assert len(list((run / "state/stages/history/qc").glob("*.json"))) == 1


def test_spatial_receipts_are_scoped_by_output_label(tmp_path, monkeypatch):
    from scoop.results import inspect_run

    run = tmp_path / "runs/test"
    (run / "config").mkdir(parents=True)
    config = run / "config/run.yaml"
    monkeypatch.setitem(sys.modules, "receipt_test", SimpleNamespace(main=lambda: None))
    for label in ("fast", "standard"):
        config.write_text(yaml.safe_dump({"spatial": {"enabled": True, "rctd": {"output_label": label}}}))
        run_command_module("spatial", "receipt_test", ["--config", str(config), "--dry-run"])
    assert sorted(path.name for path in (run / "state/stages").glob("*.json")) == ["spatial-fast.json", "spatial-standard.json"]
    assert json.loads((run / "state/stages/spatial-fast.json").read_text())["stage"] == "spatial"
    assert [stage["stage"] for stage in inspect_run(run)["stages"]] == ["spatial-fast", "spatial-standard"]


def test_atlas_approval_keeps_the_receipt_it_replaces(tmp_path):
    from scoop.cli import main

    run = tmp_path / "runs/test"
    (run / "atlas").mkdir(parents=True)
    paths = []
    for name in ("input.h5ad", "report.json", "labels.tsv", "indices.tsv"):
        (tmp_path / name).write_text(name)
        paths.append(str(tmp_path / name))
    arguments = ["atlas", "approve", "--gate", str(run / "atlas/approval.yaml"), "--run-id", "test", "--approved-by", "me",
                 "--input-h5ad", paths[0], "--annotation-report", paths[1], "--label-table", paths[2],
                 "--labeled-indices", paths[3]]
    assert main([*arguments, "--annotation-version", "v1"]) == 0
    assert main([*arguments, "--annotation-version", "v2"]) == 0
    assert json.loads((run / "state/stages/atlas-approval.json").read_text())["annotation_version"] == "v2"
    [archived] = (run / "state/stages/history/atlas-approval").glob("*.json")
    assert json.loads(archived.read_text())["annotation_version"] == "v1"
