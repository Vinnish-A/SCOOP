from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from scoop.approval import approve, verify
from scoop.atlas import report, train_scpoli


def _atlas_run(tmp_path: Path) -> tuple[Path, Path]:
    run_root = tmp_path / "runs" / "atlas-test"
    config_path = run_root / "config" / "run.yaml"
    config_path.parent.mkdir(parents=True)
    config_path.write_text(
        yaml.safe_dump(
            {
                "run": {"run_id": "atlas-test", "organism": "human", "tissue": "lung"},
                "atlas": {"graph": {"cluster_key": "leiden_scvi"}, "annotation": {"version": "atlas-v1"}},
            }
        ),
        encoding="utf-8",
    )
    obs = pd.DataFrame(
        {
            "leiden_scvi": pd.Categorical(["a", "a", "b", "b"]),
            "sample_id": ["s1", "s2", "s1", "s2"],
            "donor_id": ["d1", "d2", "d1", "d2"],
            "study_id": ["study", "study", "study", "study"],
            "condition": ["control", "case", "control", "case"],
            "scvi_curated_cell_type": ["T cell", "T cell", "B cell", "B cell"],
            "annotation_confidence": ["high", "high", "medium", "medium"],
            "annotation_status": ["labeled", "labeled", "unlabeled", "unlabeled"],
            "annotation_source": ["scvi_cluster_curated"] * 4,
            "annotation_version": ["atlas-v1"] * 4,
        },
        index=["c1", "c2", "c3", "c4"],
    )
    adata = ad.AnnData(X=sparse.csr_matrix(np.ones((4, 3))), obs=obs)
    adata.obsm["X_scVI"] = np.array([[0, 0], [0.1, 0], [2, 2], [2.1, 2]])
    adata.obsm["X_nmf_usage"] = np.array([[0.9, 0.1], [0.8, 0.2], [0.2, 0.8], [0.1, 0.9]])
    adata.uns["nmf_program_names"] = ["T_identity", "B_identity"]
    adata.obsp["scvi_connectivities"] = sparse.csr_matrix(
        np.array(
            [
                [0, 1, 0.1, 0],
                [1, 0, 0, 0.1],
                [0.1, 0, 0, 1],
                [0, 0.1, 1, 0],
            ]
        )
    )
    source = run_root / "artifacts" / "adata_scvi_annotated.h5ad"
    source.parent.mkdir(parents=True)
    adata.write_h5ad(source)
    label_dir = run_root / "atlas" / "tables"
    label_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "cluster_id": ["a", "b"],
            "scvi_curated_cell_type": ["T cell", "B cell"],
            "annotation_confidence": ["high", "medium"],
            "annotation_status": ["labeled", "unlabeled"],
            "annotation_source": ["scvi_cluster_curated"] * 2,
            "annotation_version": ["atlas-v1"] * 2,
        }
    ).to_csv(label_dir / "scvi_curated_labels.tsv", sep="\t", index=False)
    return config_path, source


def test_atlas_report_contains_distributions_neighbours_programmes_and_supervision(tmp_path: Path) -> None:
    config, source = _atlas_run(tmp_path)
    html, labels, indices = report(config_path=config, input_path=source)
    text = html.read_text(encoding="utf-8")
    assert "Subgroup distribution" in text
    assert "Graph neighbours" in text
    assert "cNMF / FastCNMF usage" in text
    assert "scPoli blocked" in text
    assert labels.is_file()
    supervised = pd.read_csv(indices, sep="\t")
    assert supervised["cell_id"].tolist() == ["c1", "c2"]


def test_approval_hashes_go_stale_and_scpoli_refuses_missing_gate(tmp_path: Path) -> None:
    config, source = _atlas_run(tmp_path)
    html, labels, indices = report(config_path=config, input_path=source)
    gate = tmp_path / "approval.yaml"
    approve(
        gate_path=gate,
        run_id="atlas-test",
        annotation_version="atlas-v1",
        approved_by="reviewer",
        input_h5ad=source,
        annotation_report=html,
        label_table=labels,
        labeled_indices=indices,
    )
    assert verify(gate)["ok"]
    html.write_text(html.read_text(encoding="utf-8") + "\nchanged", encoding="utf-8")
    result = verify(gate)
    assert not result["ok"]
    assert "approval artifact hash mismatch: annotation_report" in result["errors"]
    try:
        train_scpoli(gate_path=tmp_path / "missing.yaml")
    except PermissionError as exc:
        assert "approval gate rejected" in str(exc)
    else:
        raise AssertionError("scPoli should refuse a missing approval gate")


def test_scpoli_rejects_content_inconsistency_even_when_hashes_are_valid(tmp_path: Path) -> None:
    config, source = _atlas_run(tmp_path)
    html, labels, indices = report(config_path=config, input_path=source)
    label_table = pd.read_csv(labels, sep="\t")
    label_table.loc[label_table["cluster_id"].eq("a"), "scvi_curated_cell_type"] = "Wrong label"
    label_table.to_csv(labels, sep="\t", index=False)
    gate = tmp_path / "inconsistent-approval.yaml"
    approve(
        gate_path=gate,
        run_id="atlas-test",
        annotation_version="atlas-v1",
        approved_by="reviewer",
        input_h5ad=source,
        annotation_report=html,
        label_table=labels,
        labeled_indices=indices,
    )
    assert verify(gate)["ok"]
    try:
        train_scpoli(gate_path=gate)
    except PermissionError as exc:
        assert "label table disagrees with H5AD" in str(exc)
    else:
        raise AssertionError("scPoli should reject an internally inconsistent approval")
