from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from scsp_agent_sop.tumor_cnv import prepare_tumor_fastcnv_reference, run_tumor_fastcnv_workflow


def test_prepare_tumor_fastcnv_reference_marks_normal_pools() -> None:
    obs = pd.DataFrame(
        {
            "sample_id": ["s1"] * 4 + ["s2"] * 4,
            "cell_type_lvl1": ["T", "Myeloid", "Epithelial", "Epithelial", "B", "Endothelial", "Epithelial", "Epithelial"],
            "annotation_confidence": ["high"] * 8,
            "tumor_status": ["normal", "normal", "normal", "tumor", "normal", "normal", "mixed", "tumor"],
        },
        index=[f"c{i}" for i in range(8)],
    )
    adata = ad.AnnData(X=np.ones((8, 3)), obs=obs, var=pd.DataFrame(index=["G1", "G2", "G3"]))

    plan = prepare_tumor_fastcnv_reference(
        adata,
        major_lineage_key="cell_type_lvl1",
        sample_key="sample_id",
        parenchymal_lineages=("Epithelial",),
        normal_status_key="tumor_status",
        min_reference_cells=2,
    )

    assert plan.should_run_fastcnv
    assert plan.reference_labels == ["normal_nonparenchymal", "normal_parenchymal"]
    assert plan.n_normal_nonparenchymal == 4
    assert plan.n_normal_parenchymal == 1
    assert plan.n_candidate_cells == 3
    assert adata.obs.loc["c0", "fastcnv_reference_pool"] == "normal_nonparenchymal"
    assert adata.obs.loc["c2", "fastcnv_reference_pool"] == "normal_parenchymal"
    assert adata.obs.loc["c3", "fastcnv_reference_pool"] == "candidate_tumor_or_mixed"


def test_unresolved_low_confidence_and_filtered_cells_are_never_reference():
    import pytest
    obj = ad.AnnData(np.ones((5, 3)), obs=pd.DataFrame({
        'sample_id': ['s'] * 5,
        'cell_type_lvl1': ['Unresolved', 'T', 'T', 'T', 'Epithelial'],
        'annotation_confidence': ['high', 'low', 'high', 'high', 'high'],
        'final_use': [True, True, False, True, True],
    }, index=[str(i) for i in range(5)]))
    plan = prepare_tumor_fastcnv_reference(obj, min_reference_cells=1)
    assert plan.n_reference_cells == 1
    assert obj.obs.fastcnv_reference_pool.eq('normal_nonparenchymal').tolist() == [False, False, False, True, False]
    del obj.obs['annotation_confidence']
    with pytest.raises(ValueError, match='reviewed'):
        prepare_tumor_fastcnv_reference(obj)


def test_committed_and_legacy_statuses_respect_explicit_reference_review():
    obs = pd.DataFrame({
        'sample_id': ['s'] * 5,
        'cell_type_lvl1': ['T'] * 5,
        'annotation_status': ['accepted', 'accepted', 'labeled', 'review_required', 'unlabeled'],
        'reviewed_reference': [True, False, True, True, True],
    }, index=[f'c{i}' for i in range(5)])
    obj = ad.AnnData(np.ones((5, 3)), obs=obs)
    plan = prepare_tumor_fastcnv_reference(
        obj, reference_eligible_key='reviewed_reference', min_reference_cells=2,
    )
    assert plan.n_reference_cells == 2
    assert obj.obs.fastcnv_reference_pool.eq('normal_nonparenchymal').tolist() == [True, False, True, False, False]


def test_duplicate_symbols_cannot_silently_select_a_genomic_column(tmp_path):
    import pytest
    obj = ad.AnnData(np.array([[1, 100, 2], [3, 200, 4]]),
        var=pd.DataFrame(index=['ATRIP', 'ATRIP', 'G2']))
    with pytest.raises(ValueError, match='CNV genomic positions require unique'):
        run_tumor_fastcnv_workflow(obj, pd.DataFrame(), output_dir=tmp_path)


def _fraction_result(adata, values):
    from types import SimpleNamespace
    return SimpleNamespace(cell_metadata=pd.DataFrame({"cnv_fraction": values}, index=adata.obs_names))


def test_too_few_references_leave_calls_unassessed_unless_q75_fallback_is_opted_in():
    from scsp_agent_sop.tumor_cnv import attach_tumor_fastcnv_evidence
    obj = ad.AnnData(np.ones((6, 2)), obs=pd.DataFrame(
        {"pool": ["normal_nonparenchymal"] * 2 + ["candidate_tumor_or_mixed"] * 4}, index=[f"c{i}" for i in range(6)]))
    attach_tumor_fastcnv_evidence(obj, _fraction_result(obj, [.01, .02, .1, .2, .3, .4]), reference_key="pool")
    assert obj.obs.fastcnv_threshold_source.eq("unassessed_insufficient_reference").all()
    assert obj.obs.fastcnv_normal_threshold.isna().all()
    assert obj.obs.fastcnv_tumor_evidence.tolist() == ["normal_reference"] * 2 + ["cnv_unassessed"] * 4
    attach_tumor_fastcnv_evidence(obj, _fraction_result(obj, [.01, .02, .1, .2, .3, .4]), reference_key="pool",
                                  all_cells_q75_fallback=True)
    assert obj.obs.fastcnv_threshold_source.eq("all_cells_q75_fallback").all()
    assert obj.obs.fastcnv_tumor_evidence.tolist()[2:] == ["cnv_normal_like"] * 2 + ["cnv_tumor_like"] * 2


def test_plan_records_vocabulary_and_reviewed_cells_it_drops():
    obj = ad.AnnData(np.ones((5, 2)), obs=pd.DataFrame({
        "sample_id": ["s"] * 5,
        "cell_type_lvl1": ["Myeloid", "T_NK", "T_NK", "Oligodendrocyte", "Glioneural"],
        "annotation_confidence": ["high", "high", "high", "low", "high"],
    }, index=[f"c{i}" for i in range(5)]))
    plan = prepare_tumor_fastcnv_reference(obj, parenchymal_lineages=("Glioneural",), min_reference_cells=1)
    assert "Myeloid" in plan.nonparenchymal_lineages and "T_NK" not in plan.nonparenchymal_lineages
    assert plan.dropped_reference_by_lineage == {"T_NK": 2}
    plan = prepare_tumor_fastcnv_reference(obj, parenchymal_lineages=("Glioneural",),
                                           nonparenchymal_lineages={"Myeloid", "T_NK"}, min_reference_cells=1)
    assert plan.nonparenchymal_lineages == ["Myeloid", "T_NK"] and plan.dropped_reference_by_lineage == {}


def test_skipped_plan_is_reported_as_not_run(tmp_path, monkeypatch):
    import sys
    import yaml
    from scoop.commands import tumor_cnv as command
    run_root = tmp_path / "runs/r1"; (run_root / "config").mkdir(parents=True)
    config = run_root / "config/run.yaml"; config.write_text(yaml.safe_dump({"tumor_fastcnv": {"parenchymal_lineages": "Glioneural"}}))
    source = tmp_path / "in.h5ad"
    ad.AnnData(np.ones((3, 2)), obs=pd.DataFrame({"sample_id": ["s"] * 3, "cell_type_lvl1": ["Glioneural"] * 3,
        "annotation_confidence": ["high"] * 3}, index=["a", "b", "c"]), var=pd.DataFrame(index=["G1", "G2"])).write_h5ad(source)
    genes = tmp_path / "genes.tsv"
    pd.DataFrame({"hgnc_symbol": ["G1", "G2"], "chromosome_num": [1, 1], "start_position": [1, 10],
                  "end_position": [5, 15], "chr_arm": ["p", "p"]}).to_csv(genes, sep="\t", index=False)
    monkeypatch.setattr(sys, "argv", ["scoop tumor-cnv", "--config", str(config), "--input", str(source),
                                      "--output", str(tmp_path / "out.h5ad"), "--gene-metadata", str(genes)])
    receipt = command.main()
    assert receipt["execution"]["status"] == "not_run"
    assert receipt["reference_plan"]["n_reference_cells"] == 0 and receipt["threshold_source"] is None
    from scoop.results import result_payload
    assert result_payload(receipt, base=tmp_path, route="managed")["execution"]["status"] == "not_run"


def _measured_cnv_object():
    from scsp_agent_sop.measurement import FeatureMeasurement
    obj = ad.AnnData(np.ones((24, 3)), obs=pd.DataFrame({
        "sample_id": ["s1"] * 12 + ["s2"] * 12,
        "cell_type_lvl1": (["T"] * 11 + ["Epithelial"]) * 2,
        "annotation_confidence": ["high"] * 24,
    }, index=[f"c{i}" for i in range(24)]), var=pd.DataFrame(index=["G1", "G2", "G3"]))
    measured = np.array([[True, True, True], [True, False, True]])  # s2 did not quantify G2
    genes = pd.DataFrame({"hgnc_symbol": ["G1", "G2", "G3"], "chromosome_num": [1, 1, 2],
                          "start_position": [1, 10, 1], "end_position": [5, 15, 5], "chr_arm": ["p", "p", "q"]})
    return obj, genes, FeatureMeasurement.from_var("sample_id", ["s1", "s2"], obj.var_names, measured)


def test_declared_measurement_keeps_cnv_on_genes_every_sample_measured(tmp_path, monkeypatch):
    import scsp_agent_sop.tumor_cnv as workflow
    obj, genes, fm = _measured_cnv_object()
    seen = {}
    monkeypatch.setattr(workflow, "run_fastcnv_pooled_anndata",
                        lambda adata, gene_metadata, **kw: seen.setdefault("genes", list(gene_metadata.hgnc_symbol)))
    monkeypatch.setattr(workflow, "write_pooled_outputs", lambda **kw: None)
    monkeypatch.setattr(workflow, "attach_tumor_fastcnv_evidence", lambda *a, **kw: None)
    _, plan = run_tumor_fastcnv_workflow(obj, genes, output_dir=tmp_path, measurement=fm)
    assert seen["genes"] == ["G1", "G3"]
    assert {k: plan.measurement[k] for k in ("status", "n_samples", "n_genes_excluded", "n_positioned_genes_used")} == \
        {"status": "applied", "n_samples": 2, "n_genes_excluded": 1, "n_positioned_genes_used": 2}
    table = pd.read_csv(plan.measurement["excluded_genes_table"], sep="\t")
    assert table.to_dict("records") == [{"gene": "G2", "n_samples_unmeasured": 1, "positioned": True}]
    seen.clear()
    _, plan = run_tumor_fastcnv_workflow(obj, genes, output_dir=tmp_path)
    assert seen["genes"] == ["G1", "G2", "G3"] and plan.measurement == {"status": "not_declared"}


def test_cnv_receipt_records_declared_measurement(tmp_path, monkeypatch):
    import sys
    import yaml
    from scoop.commands import tumor_cnv as command
    from scsp_agent_sop.measurement import FORMAT
    from scsp_agent_sop.storage import init_file_registry, register_file
    obj, genes, fm = _measured_cnv_object()
    obj.obs["cell_type_lvl1"] = "Epithelial"  # no reference: the plan is skipped, the exclusion still recorded
    init_file_registry(obj, "r1")
    register_file(obj, key="feature_measurement", path=fm.save(tmp_path / "fm.npz"), category="artifacts", schema=FORMAT)
    run_root = tmp_path / "runs/r1"; (run_root / "config").mkdir(parents=True)
    config = run_root / "config/run.yaml"; config.write_text(yaml.safe_dump({"tumor_fastcnv": {"parenchymal_lineages": "Epithelial"}}))
    obj.write_h5ad(tmp_path / "in.h5ad")
    genes.to_csv(tmp_path / "genes.tsv", sep="\t", index=False)
    monkeypatch.setattr(sys, "argv", ["scoop tumor-cnv", "--config", str(config), "--input", str(tmp_path / "in.h5ad"),
                                      "--output", str(tmp_path / "out.h5ad"), "--gene-metadata", str(tmp_path / "genes.tsv")])
    receipt = command.main()
    recorded = receipt["reference_plan"]["measurement"]
    assert recorded["status"] == "applied" and recorded["n_positioned_genes_excluded"] == 1
    assert str((tmp_path / "fm.npz").resolve()) in receipt["inputs"]
    assert str(Path(recorded["excluded_genes_table"]).resolve()) in receipt["outputs"]
    assert any("positioned genes were excluded" in item for item in receipt["limitations"])
