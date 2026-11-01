"""Phase G W6: typed obs keys (R10), obs-parquet sidecars (DC1, C1-08, S032), masked COSG."""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd
import pytest
import yaml
from anndata import AnnData, read_h5ad
from scipy import sparse

from scsp_agent_sop.storage import init_file_registry, read_adata, sha256_file, write_obs_sidecar


def test_obs_sidecar_carries_declared_columns_and_overlays_only_onto_its_base(tmp_path):
    base = tmp_path / "base.h5ad"
    adata = AnnData(np.ones((3, 2)), obs=pd.DataFrame({"sample_id": ["a", "a", "b"], "final_use": [True, True, False]},
                                                      index=["c0", "c1", "c2"]))
    adata.write_h5ad(base)
    init_file_registry(adata, "r1")
    before = adata.obs.copy()
    adata.obs["final_use"] = [True, True, False]  # recomputed to its input value, still a declared output
    adata.obs["label"] = np.array(["T", None, 3], dtype=object)
    adata.obs["score"] = [0.1, 0.2, 0.3]
    receipt = write_obs_sidecar(adata, ["final_use", "label", "score"], tmp_path / "obs.parquet", base=base)
    table = pd.read_parquet(tmp_path / "obs.parquet")
    assert list(table.columns) == ["obs_name", "final_use", "label", "score"] and table.obs_name.tolist() == ["c0", "c1", "c2"]
    record = json.loads(receipt.read_text())
    assert record["sidecar"]["columns"] == ["final_use", "label", "score"] and record["sidecar"]["n_obs"] == 3
    assert record["sidecar"]["base"]["path"] == str(base.resolve()) and record["sidecar"]["base"]["sha256"] == sha256_file(base)
    assert record["file_registry"]["tables"]["obs"]["sha256"] == record["sidecar"]["sha256"]
    overlaid = read_adata(tmp_path / "obs.parquet")
    assert overlaid.obs.sample_id.tolist() == ["a", "a", "b"] and overlaid.obs.final_use.tolist() == [True, True, False]
    assert overlaid.obs.label.tolist()[0] == "T" and overlaid.uns["file_registry"]["tables"]["obs"]["columns"]
    overlaid.write_h5ad(tmp_path / "overlaid.h5ad")  # restored label columns stay writable
    with pytest.raises(ValueError, match="parquet"):
        write_obs_sidecar(adata, ["label"], tmp_path / "obs.h5ad", base=base)
    AnnData(np.zeros((3, 2)), obs=before).write_h5ad(base)  # a later run replaces the base
    with pytest.raises(ValueError, match="changed after"):
        read_adata(tmp_path / "obs.parquet")


def test_review_and_commit_read_the_sidecar_tumor_cnv_declared_not_an_older_h5ad(tmp_path, monkeypatch):
    from types import SimpleNamespace

    from scoop.commands import annotation_commit, annotation_review
    from scoop.workflow import run_command_module

    obs = pd.DataFrame({"sample_id": ["s1", "s2"] * 3, "cluster_identity": pd.Categorical(["0", "0", "0", "1", "1", "1"])},
                       index=[f"c{i}" for i in range(6)])
    run, call = _run(tmp_path, monkeypatch, {"run": {"run_id": "r1"}, "annotation": {"decision": {"min_cells_per_sample": 1}}})
    evidence = run / "artifacts" / "adata_annotation_evidence.h5ad"
    AnnData(np.ones((6, 2)), obs=obs).write_h5ad(evidence)
    AnnData(np.ones((6, 2)), obs=obs).write_h5ad(run / "artifacts" / "adata_tumor_fastcnv.h5ad")  # an earlier h5ad-mode run
    sidecar = run / "artifacts" / "tumor_fastcnv_obs.parquet"

    def tumor_cnv():
        adata = read_h5ad(evidence)
        init_file_registry(adata, "r1")
        adata.obs["fastcnv_tumor_evidence"] = np.array(["cnv_tumor_like"] * 3 + ["normal_reference"] * 3, dtype=object)
        receipt = write_obs_sidecar(adata, ["fastcnv_tumor_evidence"], sidecar, base=evidence, key="tumor_fastcnv_obs")
        return {"inputs": [str(evidence)], "outputs": [str(sidecar), str(receipt)], "output_format": "obs-parquet"}

    monkeypatch.setitem(sys.modules, "fake_tumor_cnv", SimpleNamespace(main=tumor_cnv))
    assert run_command_module("tumor-cnv", "fake_tumor_cnv", ["--config", str(run / "config" / "run.yaml")]) == 0
    assert annotation_review.default_input(run, {}) == annotation_commit.default_input(run, {}) == sidecar.resolve()
    review = call(annotation_review)
    assert review["inputs"][0] == str(sidecar.resolve())
    out = run / "annotation" / "decisions"
    assert json.loads((out / "annotation_evidence_bundle.json").read_text())["clusters"][0]["cnv_refs"] == ["fastcnv_tumor_evidence"]
    template = json.loads((out / "annotation_decision_template.json").read_text())["decisions"]
    decisions = out / "decisions.json"
    decisions.write_text(json.dumps([{**template[0], "parent_label": "Epithelial", "canonical_label": "Epithelial",
                                      "final_label": "Epithelial", "confidence": "low", "review_required": False,
                                      "reason": "synthetic"}]))
    commit = call(annotation_commit, "--decisions", str(decisions), "--output-format", "obs-parquet")
    assert commit["inputs"][0] == str(sidecar.resolve())
    committed = read_adata(run / "artifacts" / "annotation_committed_obs.parquet")  # sidecar on a sidecar
    assert committed.obs.fastcnv_tumor_evidence.tolist()[:3] == ["cnv_tumor_like"] * 3
    assert committed.obs.annotation_reviewed.tolist() == [True] * 3 + [False] * 3
    call(annotation_commit, "--decisions", str(decisions))
    assert "fastcnv_tumor_evidence" in read_h5ad(run / "artifacts" / "adata_annotation_committed.h5ad").obs


def _run(tmp_path, monkeypatch, cfg):
    run = tmp_path / "runs" / "test"
    (run / "config").mkdir(parents=True)
    (run / "artifacts").mkdir()
    config = run / "config" / "run.yaml"
    config.write_text(yaml.safe_dump(cfg))
    monkeypatch.chdir(tmp_path)

    def call(module, *args):
        monkeypatch.setattr(sys, "argv", ["scoop", "--config", str(config), *args])
        return module.main()
    return run, call


def test_review_and_commit_take_typed_keys_a_cluster_table_and_write_an_obs_sidecar(tmp_path, monkeypatch):
    from scoop.commands import annotation_commit, annotation_review
    from scsp_agent_sop.storage import register_file

    obs = pd.DataFrame({"lib": ["s1", "s1", "s2", "s1", "s2", "s3"]}, index=[f"c{i}" for i in range(6)])
    adata = AnnData(np.ones((6, 2)), obs=obs)
    pd.DataFrame({"group": ["0", "0"], "names": ["CD3E", "CD2"], "scores": [9.0, 5.0]}).to_parquet(tmp_path / "deg.parquet")
    register_file(adata, key="cluster_markers", path=tmp_path / "deg.parquet", schema="cluster_markers.v1",
                  extra={"label_key": "cluster_r1"})
    run, call = _run(tmp_path, monkeypatch, {"run": {"run_id": "run1", "obs_keys": {"sample": "lib"}},
                                             "annotation": {"decision": {"cluster_key": "cluster_r1", "min_cells_per_sample": 1,
                                                                         "min_cells_high_confidence": 3,
                                                                         "min_samples_high_confidence": 2}}})
    adata.write_h5ad(run / "artifacts" / "adata_annotation_evidence.h5ad")  # no cluster column (C1-08)
    clusters = tmp_path / "clusters.tsv"
    pd.DataFrame({"obs_name": [f"c{i}" for i in range(6)][::-1], "cluster_r1": ["1", "1", "1", "0", "0", "0"]}).to_csv(
        clusters, sep="\t", index=False)
    result = call(annotation_review, "--cluster-table", str(clusters))
    assert result["obs_keys"]["sample"] == {"key": "lib", "source": "run.obs_keys.sample"}
    out = run / "annotation" / "decisions"
    bundle = json.loads((out / "annotation_evidence_bundle.json").read_text())
    assert bundle["composition_keys"]["sample_key"] == "lib"
    template = json.loads((out / "annotation_decision_template.json").read_text())["decisions"]
    decision = {**template[0], "parent_label": "Lymphoid", "canonical_label": "T", "final_label": "T_cell",
                "confidence": "medium", "positive_markers": ["CD3E"], "review_required": False, "reason": "synthetic"}
    decisions = out / "decisions.json"
    decisions.write_text(json.dumps([decision]))
    receipt = call(annotation_commit, "--decisions", str(decisions), "--cluster-table", str(clusters),
                   "--output-format", "obs-parquet")
    sidecar = run / "artifacts" / "annotation_committed_obs.parquet"
    assert not (run / "artifacts" / "adata_annotation_committed.h5ad").exists()
    assert str(sidecar.resolve()) in receipt["outputs"] and str(clusters.resolve()) in receipt["inputs"]
    table = pd.read_parquet(sidecar).set_index("obs_name")
    assert "cluster_r1" not in table and "lib" not in table
    # the table lists obs in reverse, so cluster 0 is c0-c2
    assert table.loc[["c0", "c1", "c2"], "annotation_status"].tolist() == ["accepted"] * 3
    assert table.loc["c3", "annotation_status"] == "unassigned" and not table.loc["c3", "annotation_reviewed"]
    registry = json.loads(sidecar.with_suffix(".json").read_text())["file_registry"]
    assert {"annotation_decisions", "annotation_committed_obs"} <= set(registry["tables"])

    with pytest.raises(ValueError, match="no 'cluster_r1' for 1 of 6"):
        partial = tmp_path / "partial.tsv"
        pd.read_csv(clusters, sep="\t").iloc[1:].to_csv(partial, sep="\t", index=False)
        call(annotation_commit, "--decisions", str(decisions), "--cluster-table", str(partial))


def test_tumor_cnv_obs_parquet_writes_no_h5ad(tmp_path, monkeypatch):
    from scoop.commands import tumor_cnv as command
    import scsp_agent_sop.tumor_cnv as workflow

    obs = pd.DataFrame({"sample_id": ["s1"] * 12 + ["s2"] * 12, "cell_type_lvl1": (["T"] * 11 + ["Epithelial"]) * 2,
                        "annotation_confidence": ["high"] * 24}, index=[f"c{i}" for i in range(24)])
    source = tmp_path / "in.h5ad"
    AnnData(np.ones((24, 2)), obs=obs, var=pd.DataFrame(index=["G1", "G2"])).write_h5ad(source)
    genes = tmp_path / "genes.tsv"
    pd.DataFrame({"hgnc_symbol": ["G1", "G2"], "chromosome_num": [1, 1], "start_position": [1, 10],
                  "end_position": [5, 15], "chr_arm": ["p", "p"]}).to_csv(genes, sep="\t", index=False)

    def fake_cnv(adata, gene_metadata, **kw):
        adata.obs["fastcnv_cnv_fraction"] = 0.1
        adata.obs["fastcnv_threshold_source"] = "reference_q99_mad"
        return {}
    monkeypatch.setattr(workflow, "run_fastcnv_pooled_anndata", fake_cnv)
    monkeypatch.setattr(workflow, "write_pooled_outputs", lambda **kw: None)
    monkeypatch.setattr(workflow, "attach_tumor_fastcnv_evidence", lambda *a, **kw: None)
    run, call = _run(tmp_path, monkeypatch, {"tumor_fastcnv": {"parenchymal_lineages": "Epithelial"}})
    receipt = call(command, "--input", str(source), "--gene-metadata", str(genes), "--output-format", "obs-parquet")
    sidecar = run / "artifacts" / "tumor_fastcnv_obs.parquet"
    assert receipt["output_format"] == "obs-parquet" and str(sidecar.resolve()) in receipt["outputs"]
    assert not (run / "artifacts" / "adata_tumor_fastcnv.h5ad").exists()
    table = pd.read_parquet(sidecar)
    assert {"fastcnv_cnv_fraction", "fastcnv_threshold_source"} <= set(table.columns)
    assert not {"sample_id", "cell_type_lvl1", "annotation_confidence"} & set(table.columns)


def test_qc_obs_parquet_needs_an_explicit_output_and_writes_qc_columns(tmp_path, monkeypatch):
    from scoop.commands import qc

    a = AnnData(sparse.csr_matrix([[6, 7, 99], [1, 0, 0], [22, 11, 1], [30, 2, 0]]),
                obs=pd.DataFrame({"sample_id": ["A"] * 4}, index=[f"c{i}" for i in range(4)]))
    a.var_names = ["CD1C", "FCER1A", "NegativeControlProbe_1"]
    a.var["feature_types"] = ["Gene Expression", "Gene Expression", "Negative Control Probe"]
    a.layers["counts"] = a.X.copy()
    a.obsm["spatial"] = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0]])
    run, call = _run(tmp_path, monkeypatch, {"run": {"assay": "xenium"}, "qc": {"spatial": {"min_counts": 10, "min_genes": 2}}})
    a.write_h5ad(run / "artifacts" / "validated_adata.h5ad")
    with pytest.raises(SystemExit):
        call(qc, "--output-format", "obs-parquet")
    call(qc, "--output-format", "obs-parquet", "--output", str(run / "artifacts" / "qc_obs.parquet"))
    assert not (run / "artifacts" / "adata_qc.h5ad").exists()
    table = pd.read_parquet(run / "artifacts" / "qc_obs.parquet")
    assert table.final_use.tolist() == [True, False, True, True] and "sample_id" not in table
    assert "qc_thresholds_by_sample" in json.loads((run / "artifacts" / "qc_obs.json").read_text())["file_registry"]["tables"]


def _programme_object(tmp_path):
    rng = np.random.default_rng(0)
    counts = rng.poisson(2.0, size=(60, 40)).astype(np.float32)
    counts[:30, :10] += rng.poisson(6.0, size=(30, 10))
    counts[30:, 10:20] += rng.poisson(6.0, size=(30, 10))
    counts[::2] *= 2  # a technical batch effect on half of the cells
    adata = AnnData(sparse.csr_matrix(counts), obs=pd.DataFrame({"sample_id": ["s1", "s2", "s3"] * 20,
                                                                  "lane": ["L1", "L2"] * 30},
                                                                 index=[f"c{i}" for i in range(60)]),
                    var=pd.DataFrame(index=[f"G{i}" for i in range(40)]))
    adata.layers["counts"] = adata.X.copy()
    adata.layers["log1p_norm"] = sparse.csr_matrix(np.log1p(counts))
    return adata


def test_programs_pass_the_typed_batch_to_fastcnmf_preprocess(tmp_path, monkeypatch):
    from scoop.commands import programs

    nmf = {"k_grid": [2], "seeds": [0, 1], "max_iter": 50, "workers": 1, "n_top_hvg": 20}
    run, call = _run(tmp_path, monkeypatch, {"run": {"obs_keys": {"batch": "lane"}}, "programs": nmf})
    _programme_object(tmp_path).write_h5ad(run / "artifacts" / "adata_core.h5ad")
    receipt = call(programs)
    assert receipt["batch_handling"] == {"method": "none", "key": "lane", "source": "run.obs_keys.batch"}
    assert any("did not model the declared batch 'lane'" in item for item in receipt["limitations"])

    run.joinpath("config", "run.yaml").write_text(yaml.safe_dump({"run": {"obs_keys": {"batch": "lane"}},
                                                                  "programs": {**nmf, "batch_correction": "harmony_moe"}}))
    receipt = call(programs)
    manifest = json.loads((run / "programs" / "fastcnmf_preprocess" / "fastcnmf_preprocess_manifest.json").read_text())
    assert manifest["sample_key"] == "lane" and manifest["n_vars"] == 20
    assert receipt["programme_method"]["batch_key"] == "lane"
    assert receipt["programme_method"]["expression"] == "fastcnmf_harmony_moe_corrected_hvg"
    out = read_h5ad(run / "artifacts" / "adata_programs.h5ad")
    assert out.obsm["X_nmf_usage"].shape == (60, 2) and int(out.var.nmf_discovery_gene.sum()) == 20
    assert manifest["corrected_h5ad"] in receipt["outputs"]

    run.joinpath("config", "run.yaml").write_text(yaml.safe_dump({"programs": {**nmf, "batch_correction": "harmony_moe"},
                                                                  "keys": {"batch_candidates": ["chemistry"]}}))
    with pytest.raises(KeyError, match=r"batch=None"):
        call(programs)


def test_cosg_scores_only_genes_every_sample_measured(tmp_path, monkeypatch):
    from scoop.commands import annotation_evidence
    from scsp_agent_sop.annotation import run_markers_fastde
    from scsp_agent_sop.measurement import FORMAT, FeatureMeasurement
    from scsp_agent_sop.storage import register_file

    matrix = np.log1p(np.array([[6, 5, 0, 0], [6, 4, 0, 1], [0, 0, 5, 6], [1, 0, 6, 5]] * 2, dtype=float))
    adata = AnnData(matrix, obs=pd.DataFrame({"cluster_identity": pd.Categorical(["0", "0", "1", "1"] * 2),
                                              "sample_id": ["s1"] * 4 + ["s2"] * 4}, index=[f"c{i}" for i in range(8)]),
                    var=pd.DataFrame(index=["g0", "g1", "g2", "g3"]))
    adata.layers["log1p_norm"] = matrix.copy()
    fm = FeatureMeasurement.from_var("sample_id", ["s1", "s2"], adata.var_names,
                                     [[True] * 4, [True, False, True, True]])  # s2 did not quantify g1
    scores = run_markers_fastde(adata, layer="log1p_norm", n_genes=4, measurement=fm)
    assert "g1" not in set(scores.astype(str).to_numpy().ravel()) and adata.n_vars == 4
    assert {"g0", "g2", "g3"} & set(scores.astype(str).to_numpy().ravel())

    init_file_registry(adata, "r1")
    register_file(adata, key="feature_measurement", path=fm.save(tmp_path / "fm.npz"), category="artifacts", schema=FORMAT)
    run, call = _run(tmp_path, monkeypatch, {"run": {"run_id": "r1"}})
    adata.write_h5ad(run / "artifacts" / "adata_programs.h5ad")
    receipt = call(annotation_evidence, "--use-fastde-cosg")
    specificity = receipt["measurement"]["specificity"]
    assert specificity["gene_axis"] == "genes_measured_in_every_sample" and specificity["n_genes_excluded"] == 1
    assert pd.read_csv(specificity["excluded_genes_table"], sep="\t").gene.tolist() == ["g1"]
