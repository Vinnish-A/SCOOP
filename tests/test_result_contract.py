import json
import sys
from types import SimpleNamespace

import pytest
import yaml

from scoop.results import inspect_run, register_result, result_payload
from scoop.workflow import run_command_module


def _manifest(tmp_path):
    (tmp_path / "counts").write_bytes(b"measured")
    (tmp_path / "weights").write_bytes(b"computed")
    path = tmp_path / "result.json"
    path.write_text(json.dumps({"stage": "spatial", "inputs": ["counts"], "outputs": ["weights"],
                              "validation": {"status": "limited"},
                              "coverage": [{"scope": "spatial", "unit": "spots", "total": 10, "processed": 8}]}))
    return path


def test_registered_result_detects_changed_input_and_preserves_previous_receipt(tmp_path):
    manifest = _manifest(tmp_path)
    register_result(tmp_path, manifest)
    report = inspect_run(tmp_path, verify=True)
    assert report["stages"][0]["integrity"] == "verified"
    assert report["stages"][0]["interpretation"]["status"] == "not_assessed"
    (tmp_path / "counts").write_bytes(b"replaced")
    assert inspect_run(tmp_path, verify=True)["stale_or_missing_stages"] == ["spatial"]
    register_result(tmp_path, manifest)
    assert len(list((tmp_path / "state/stages/history/spatial").glob("*.json"))) == 1


def test_registered_not_run_result_is_skipped_not_executed(tmp_path):
    manifest = _manifest(tmp_path)
    payload = json.loads(manifest.read_text())
    manifest.write_text(json.dumps({**payload, "outputs": [], "execution": {"status": "not_run", "reason": "no reference"}}))
    receipt = json.loads(register_result(tmp_path, manifest).read_text())
    assert receipt["status"] == "skipped" and receipt["execution"]["reason"] == "no reference"


def test_passed_validation_requires_bound_inputs_and_actual_checks(tmp_path):
    with pytest.raises(ValueError, match="passed checks"):
        result_payload({"validation": {"status": "passed"}}, base=tmp_path, route="managed")
    with pytest.raises(FileNotFoundError):
        result_payload({"outputs": ["phantom"]}, base=tmp_path, route="managed")


def test_legacy_committed_does_not_imply_validation_or_interpretation(tmp_path):
    stage = tmp_path / "state/stages/core.json"
    stage.parent.mkdir(parents=True)
    stage.write_text(json.dumps({"status": "committed", "inputs": {"old": "counts"}}))
    report = inspect_run(tmp_path, verify=True)["stages"][0]
    assert report["execution"]["status"] == "completed"
    assert report["validation"]["status"] == "not_assessed"
    assert report["interpretation"]["status"] == "not_assessed"
    assert report["integrity"] == "not_fully_bound"


def test_zero_exit_and_stale_directory_are_not_validation(tmp_path, monkeypatch):
    run = tmp_path / "runs/test"
    (run / "config").mkdir(parents=True)
    config = run / "config/run.yaml"
    config.write_text(yaml.safe_dump({"run": {"run_id": "test"}, "spatial": {"enabled": True}}))
    (run / "spatial").mkdir()
    (run / "spatial/old_weights").write_bytes(b"old")
    monkeypatch.setitem(sys.modules, "receipt_contract_test", SimpleNamespace(main=lambda: None))
    assert run_command_module("spatial", "receipt_contract_test", ["--config", str(config)]) == 0
    receipt = json.loads((run / "state/stages/spatial.json").read_text())
    assert receipt["status"] == "executed"
    assert receipt["validation"]["status"] == "not_assessed"
    assert receipt["outputs"] == []
    assert inspect_run(run, verify=True)["stages"][0]["integrity"] == "not_fully_bound"


def test_help_does_not_overwrite_analysis_receipt(tmp_path, monkeypatch):
    manifest = _manifest(tmp_path)
    receipt = register_result(tmp_path, manifest)
    before = receipt.read_bytes()
    def help_main():
        raise SystemExit(0)
    monkeypatch.setitem(sys.modules, "receipt_contract_help", SimpleNamespace(main=help_main))
    assert run_command_module("spatial", "receipt_contract_help", ["--config", str(tmp_path / "missing.yaml"), "--help"]) == 0
    assert receipt.read_bytes() == before


def test_main_help_lists_semantic_commands(capsys):
    from scoop.cli import main
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0
    output = capsys.readouterr().out
    for action in ("spatial", "annotate", "programs", "deliver", "register-result", "status"):
        assert action in output
