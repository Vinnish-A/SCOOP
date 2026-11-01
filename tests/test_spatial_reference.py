import json

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse

from scoop.io import sha256_file
from scsp_agent_sop import rctd_adapter as adapter
from scsp_agent_sop.measurement import FeatureMeasurement
from scsp_agent_sop.reference_profiles import build_reference
from scsp_agent_sop.spatial_identifiability import run_identifiability

PROFILES = np.array([[.6, .1, .1], [.2, .7, .1], [.1, .1, .7], [.1, .1, .1]])  # genes g1..g3, heldout h


def _fake_fit(calls):
    def fake_fit(spatial, status, reference, options, mode):
        calls.append(True)
        y = spatial.X.toarray() if sparse.issparse(spatial.X) else np.asarray(spatial.X)
        genes = [spatial.var_names.get_loc(g) for g in reference.gene_names]
        # Least squares on the fit genes stands in for RCTD.
        raw = np.clip(np.linalg.lstsq(reference.profiles, (y[:, genes] / y[:, genes].sum(1, keepdims=True)).T, rcond=None)[0].T, 0, None)
        return {"raw": raw, "fitted_rows": np.arange(spatial.n_obs), "profiles": reference.profiles,
                "fit_genes": np.array(reference.gene_names), "types": np.array(reference.cell_type_names)}, {
                "method": "test fixture", "version": "0.4.0", "mode": mode, "sigma": 100}
    return fake_fit


def _package(tmp_path, monkeypatch, *, v2=True, pending=False, types=("A", "B", "C"), parents=("P", "P", "C")):
    rng = np.random.default_rng(0)
    truth = rng.dirichlet(np.ones(3), size=40)
    counts = rng.poisson(1000 * truth @ PROFILES.T)
    obj = AnnData(sparse.csr_matrix(counts), obs=pd.DataFrame(index=[f"s{i}" for i in range(40)]),
                  var=pd.DataFrame(index=["g1", "g2", "g3", "h"]))
    path = tmp_path / "spatial.h5ad"; obj.write_h5ad(path)
    profiles = tmp_path / "profiles.npz"
    np.savez(profiles, profiles=PROFILES, gene_names=np.array(["g1", "g2", "g3", "h"]), cell_type_names=np.array(types))
    manifest = {"unit": "mean_counts_per_full_native_exposure", "native_full_exposure_preserved": True,
                "profiles_sha256": sha256_file(profiles)}
    if v2:
        manifest.update({"schema_version": "scoop.reference-profiles.v2", "components": [
            {"name": t, "parent": p, "compartment": "tissue", "training_cells": 50,
             "review_status": "evidence_pending" if pending and t == "C" else "reviewed"} for t, p in zip(types, parents)],
            "heldout_validation_genes": {"path": "heldout.tsv"}})
        pd.DataFrame({"gene": ["h"], "selection": ["random_expressed_nonDE"], "source_component": [""]}).to_csv(
            tmp_path / "heldout.tsv", sep="\t", index=False)
    else:
        manifest.update({"schema_version": "scoop.reference-profiles.v1", "labels_reviewed": not pending})
    (tmp_path / "reference.json").write_text(json.dumps(manifest))
    lookup = tmp_path / "q.npz"; lookup.write_bytes(b"test fixture lookup")
    cfg = {"spatial": {"rctd": {"spatial_counts_layer": "X", "reference_profiles": str(profiles),
                                "reference_manifest": str(tmp_path / "reference.json"), "q_matrices": str(lookup),
                                "exclude_genes": ["h"]}}}
    calls = []
    monkeypatch.setattr(adapter, "_fit_payload", _fake_fit(calls))
    return path, cfg, calls


def test_v1_manifest_is_still_accepted_as_all_reviewed(tmp_path, monkeypatch):
    path, cfg, _ = _package(tmp_path, monkeypatch, v2=False)
    info = adapter.run_rctd_api(path, None, tmp_path / "out/r.h5ad", cfg, "full")
    summary = pd.read_parquet(tmp_path / "out/full_component_summary.parquet")
    assert summary.component.tolist() == ["A", "B", "C"] and set(summary.review_status) == {"reviewed"}
    assert info["biological_interpretation_allowed"] and info["run_purpose"] == "production"


def test_unreviewed_v1_and_evidence_pending_v2_are_refused_in_production_only(tmp_path, monkeypatch):
    (tmp_path / "v1").mkdir()
    path, cfg, calls = _package(tmp_path / "v1", monkeypatch, v2=False, pending=True)
    with pytest.raises(ValueError, match="evidence_pending"):
        adapter.run_rctd_api(path, None, tmp_path / "v1/out/r.h5ad", cfg, "full")
    (tmp_path / "v2").mkdir()
    path, cfg, calls = _package(tmp_path / "v2", monkeypatch, pending=True)
    with pytest.raises(ValueError, match=r"evidence_pending: \['C'\]"):
        adapter.run_rctd_api(path, None, tmp_path / "v2/out/r.h5ad", cfg, "full")
    assert not calls
    cfg["spatial"]["rctd"]["run_purpose"] = "pilot"
    pilot = adapter.run_rctd_api(path, None, tmp_path / "v2/out/r.h5ad", cfg, "full")
    assert pilot["biological_interpretation_allowed"] is False and pilot["evidence_pending_components"] == ["C"]
    cfg["spatial"]["rctd"]["run_purpose"] = "sensitivity"
    sensitivity = adapter.run_rctd_api(path, None, tmp_path / "v2/out/r.h5ad", cfg, "full")
    # The purpose labels the fit; it is not part of the computation.
    assert sensitivity["fit_identity"] == pilot["fit_identity"] and sensitivity["computation_reused"] and len(calls) == 1
    assert sensitivity["biological_interpretation_allowed"] is False


def test_parent_weights_are_exact_sums_and_summary_has_one_row_per_component(tmp_path, monkeypatch):
    path, cfg, _ = _package(tmp_path, monkeypatch)
    table = tmp_path / "ident.tsv"
    pd.DataFrame({"component": ["A", "B"], "r_true_vs_est": [.9, .3],
                  "identifiability_verdict": ["identifiable", "not_identifiable_fine"]}).to_csv(table, sep="\t", index=False)
    cfg["spatial"]["rctd"].update({"merge_level": "parent", "identifiability_table": str(table)})
    info = adapter.run_rctd_api(path, None, tmp_path / "out/r.h5ad", cfg, "full")
    fine = pd.read_parquet(tmp_path / "out/full_weights_visualization.parquet").set_index("spatial_unit_id")
    parent = pd.read_parquet(tmp_path / "out/full_weights_parent.parquet").set_index("spatial_unit_id")
    assert parent.columns.tolist() == ["P", "C"]
    np.testing.assert_array_equal(parent.P.to_numpy(), (fine.A + fine.B).to_numpy())
    np.testing.assert_array_equal(parent.C.to_numpy(), fine.C.to_numpy())
    summary = pd.read_parquet(tmp_path / "out/full_component_summary.parquet")
    assert len(summary) == 3 and summary.component.is_unique
    assert summary.identifiability_verdict.tolist() == ["identifiable", "not_identifiable_fine", "not_assessed"]
    np.testing.assert_allclose(summary.mean_weight, fine.mean())
    assert info["merge_level"] == "parent"


def test_heldout_genes_are_validated_against_a_shuffled_null(tmp_path, monkeypatch):
    path, cfg, _ = _package(tmp_path, monkeypatch)
    info = adapter.run_rctd_api(path, None, tmp_path / "out/r.h5ad", cfg, "full")
    table = pd.read_parquet(tmp_path / "out/full_heldout_gene_validation.parquet")
    assert table.gene.tolist() == ["h"] and table.selection.tolist() == ["random_expressed_nonDE"]
    summary = info["heldout_gene_validation"]
    assert summary["status"] == "assessed" and summary["n_genes_scored"] == 1
    assert "by_selection" in summary and np.isfinite(summary["median_shuffled_null"])


def _atlas(tmp_path, n_per_type=(60, 60, 60), donors=3, extra_obs=None):
    rng = np.random.default_rng(1)
    profiles = np.array([[.5, .1, .1], [.2, .6, .1], [.1, .1, .6], [.1, .1, .1], [.1, .1, .1]])
    labels = np.repeat(["T0", "T1", "T2"], n_per_type)
    x = np.vstack([rng.poisson(rng.integers(200, 400) * profiles[:, int(t[1])]) for t in labels])
    obs = pd.DataFrame({"cell_type": labels, "annotation_status": "accepted", "annotation_reviewed": True,
                        "donor": [f"d{i % donors}" for i in range(len(labels))],
                        "sample": [f"s{i % 2}" for i in range(len(labels))]}, index=[f"c{i}" for i in range(len(labels))])
    for key, values in (extra_obs or {}).items():
        obs[key] = values
    obj = AnnData(sparse.csr_matrix(x), obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(5)]))
    path = tmp_path / "atlas.h5ad"; obj.write_h5ad(path)
    return path, obj


def test_builder_equals_rctd_reference_profiles(tmp_path):
    import rctd

    path, obj = _atlas(tmp_path)
    manifest = build_reference(path, tmp_path / "ref", label_key="cell_type", donor_key="donor", holdout_fraction=0,
                               heldout_genes_per_component=0, heldout_random_genes=0, min_umi=100, workers=1, chunk_rows=7)
    with np.load(tmp_path / "ref/reference_profiles.npz") as z:
        ours = pd.DataFrame(z["profiles"], index=z["gene_names"], columns=z["cell_type_names"])
    expected = rctd.Reference(obj, cell_type_col="cell_type", cell_min=1, n_max_cells=10_000, min_UMI=100)
    expected = pd.DataFrame(expected.profiles, index=expected.gene_names, columns=expected.cell_type_names)
    np.testing.assert_allclose(ours, expected.loc[ours.index, ours.columns], rtol=0, atol=1e-15)
    assert manifest["schema_version"] == "scoop.reference-profiles.v2"
    assert [c["training_cells"] for c in manifest["components"]] == [60, 60, 60]
    assert {c["holdout_status"] for c in manifest["components"]} == {"not_requested"}


def test_build_reference_cli_with_workers_matches_inline_build(tmp_path):
    from scoop.cli import main

    path, _ = _atlas(tmp_path)
    build_reference(path, tmp_path / "inline", label_key="cell_type", donor_key="donor", holdout_fraction=0.2,
                    holdout_min_cells=5, heldout_random_genes=0)
    assert main(["spatial", "build-reference", "--h5ad", str(path), "--output", str(tmp_path / "cli"),
                 "--label-key", "cell_type", "--donor-key", "donor", "--holdout-fraction", "0.2",
                 "--holdout-min-cells", "5", "--heldout-random-genes", "0", "--workers", "2", "--chunk-rows", "11"]) == 0
    with np.load(tmp_path / "inline/reference_profiles.npz") as a, np.load(tmp_path / "cli/reference_profiles.npz") as b:
        np.testing.assert_allclose(a["profiles"], b["profiles"], rtol=1e-14)
    assert json.loads((tmp_path / "cli/reference_manifest.json").read_text())["holdout_rule"]["n_cells"] == 36


def test_measured_denominators_and_technical_cells_never_admitted(tmp_path):
    n = 180
    status = np.where(np.arange(n) % 9 == 0, "excluded_technical", "accepted")
    status[np.arange(n) % 9 == 1] = "review_required"
    path, obj = _atlas(tmp_path, extra_obs={"annotation_status": status})
    # Sample s1 did not measure g0: its zeros there are absence of measurement, not expression.
    x = obj.X.toarray()
    x[(obj.obs["sample"] == "s1").to_numpy(), 0] = 0
    obj.X = sparse.csr_matrix(x)
    obj.write_h5ad(path)
    measured = np.ones((2, 5), bool); measured[1, 0] = False
    measurement = FeatureMeasurement.from_var("sample", ["s0", "s1"], obj.var_names, measured)
    manifest = build_reference(path, tmp_path / "ref", label_key="cell_type", donor_key="donor", measurement=measurement,
                               admit_status=("accepted", "excluded_technical"), holdout_fraction=0,
                               heldout_genes_per_component=0, heldout_random_genes=0, min_measured_cells=1)
    x = obj.X.toarray().astype(float)
    keep = (status == "accepted") & (x.sum(1) >= 100)
    rate = x / x.sum(1, keepdims=True)
    with np.load(tmp_path / "ref/reference_profiles.npz") as z:
        ours = pd.DataFrame(z["profiles"], index=z["gene_names"], columns=z["cell_type_names"])
    for t in ours.columns:
        cells = keep & (obj.obs.cell_type == t).to_numpy()
        measuring = cells & (obj.obs["sample"] == "s0").to_numpy()
        np.testing.assert_allclose(ours.loc["g0", t], rate[measuring, 0].sum() / measuring.sum(), rtol=1e-14)
        np.testing.assert_allclose(ours.loc["g1", t], rate[cells, 1].mean(), rtol=1e-14)
        component = next(c for c in manifest["components"] if c["name"] == t)
        assert component["measured_fraction_min"] == pytest.approx(measuring.sum() / cells.sum())
    assert manifest["admission"]["n_training_cells"] == keep.sum()
    coverage = pd.read_csv(tmp_path / "ref/reference_coverage.tsv", sep="\t").set_index("label")
    assert coverage.n_excluded_technical.sum() == (status == "excluded_technical").sum()
    assert coverage.n_training.sum() == keep.sum()


def test_holdout_rule_flags_tiny_and_single_donor_components(tmp_path):
    donor = np.r_[np.full(20, "d0"), np.full(100, "d0"), [f"d{i % 3}" for i in range(100)]]
    path, _ = _atlas(tmp_path, n_per_type=(20, 100, 100), extra_obs={"donor": donor})
    manifest = build_reference(path, tmp_path / "ref", label_key="cell_type", donor_key="donor", holdout_fraction=0.2,
                               holdout_max_cells=200, holdout_min_cells=10, min_donors=2,
                               heldout_genes_per_component=0, heldout_random_genes=0)
    components = {c["name"]: c for c in manifest["components"]}
    assert components["T0"]["holdout_status"] == "below_min_cells" and components["T0"]["holdout_cells"] == 0
    assert components["T1"]["holdout_status"] == "below_min_donors" and components["T1"]["holdout_cells"] == 20
    assert components["T2"]["holdout_status"] == "sufficient" and components["T2"]["holdout_donors"] == 3
    assert components["T2"]["training_cells"] == 80
    assert any("T0: holdout below_min_cells" in gap for gap in manifest["gaps"])
    assert manifest["holdout_rule"]["min_cells"] == 10


def test_identifiability_harness_runs_on_a_built_reference(tmp_path, monkeypatch):
    path, _ = _atlas(tmp_path, n_per_type=(100, 100, 100))
    build_reference(path, tmp_path / "ref", label_key="cell_type", donor_key="donor", holdout_fraction=0.2,
                    holdout_min_cells=10, heldout_genes_per_component=0, heldout_random_genes=1)
    lookup = tmp_path / "q.npz"; lookup.write_bytes(b"test fixture lookup")
    calls = []
    monkeypatch.setattr(adapter, "_fit_payload", _fake_fit(calls))
    summary = run_identifiability(tmp_path / "ref", tmp_path / "ident", rctd_options={"q_matrices": str(lookup)},
                                  spots_per_component=10, depth_median=800, depth_sdlog=0.1)
    per = pd.read_parquet(tmp_path / "ident/component_identifiability.parquet")
    pairs = pd.read_parquet(tmp_path / "ident/pair_identifiability.parquet")
    assert per.component.tolist() == ["T0", "T1", "T2"] and len(pairs) == 3 and len(calls) == 1
    assert set(per.identifiability_verdict) <= {"identifiable", "weakly_identifiable", "not_identifiable_fine"}
    assert summary["n_fitted"] == 30 and per.dominant_recall.notna().all()
