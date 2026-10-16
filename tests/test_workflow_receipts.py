from __future__ import annotations

import json
import sys
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
