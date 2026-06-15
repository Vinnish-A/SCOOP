from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from fastcopykat import FastCopyKatConfig, run_fastcopykat
from fastcopykat.io import write_copykat_outputs


def _synthetic_copykat_input() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(7)
    n_genes = 240
    genes = [f"gene{i:03d}" for i in range(n_genes)]
    cells = [f"dip{i}" for i in range(8)] + [f"tumor{i}" for i in range(8)]
    counts = rng.poisson(5, size=(n_genes, len(cells))).astype(float)
    counts[:80, 8:] *= 3.0
    counts[160:220, 8:] *= 0.25
    counts = pd.DataFrame(np.round(counts), index=genes, columns=cells)
    chrom = np.repeat(["1", "2", "3"], n_genes // 3)
    annotation = pd.DataFrame(
        {
            "gene": genes,
            "chrom": chrom,
            "start": np.tile(np.arange(1, 81) * 1000, 3),
            "end": np.tile(np.arange(1, 81) * 1000 + 99, 3),
        }
    )
    bins = pd.DataFrame(
        {
            "chrom": np.repeat(["1", "2", "3"], 4),
            "chrompos": np.tile([20_000, 40_000, 60_000, 80_000], 3),
        }
    )
    return counts, annotation, bins


def test_fastcopykat_keeps_copykat_main_output_contract(tmp_path: Path) -> None:
    counts, annotation, bins = _synthetic_copykat_input()
    result = run_fastcopykat(
        counts,
        annotation,
        bins=bins,
        normal_cell_names=[f"dip{i}" for i in range(8)],
        sample_name="tiny",
        config=FastCopyKatConfig(
            min_gene_per_cell=10,
            min_gene_per_chromosome=1,
            low_detection_rate=0.0,
            upper_detection_rate=0.0,
            min_cluster_cells=2,
        ),
    )

    assert list(result.prediction.columns) == ["cell.names", "copykat.pred"]
    assert set(result.prediction["copykat.pred"]).issubset({"diploid", "aneuploid"})
    assert result.cna.columns[:3].tolist() == ["chrom", "chrompos", "abspos"]
    assert result.cna.shape[0] == bins.shape[0]
    assert (result.prediction.set_index("cell.names").loc[[f"dip{i}" for i in range(8)], "copykat.pred"] == "diploid").all()
    assert (result.prediction.set_index("cell.names").loc[[f"tumor{i}" for i in range(8)], "copykat.pred"] == "aneuploid").sum() >= 6

    outputs = write_copykat_outputs(
        prediction=result.prediction,
        cna=result.cna,
        output_dir=tmp_path,
        sample_name="tiny",
        manifest=result.manifest,
    )
    assert Path(outputs["prediction"]).name == "tiny_copykat_prediction.txt"
    assert Path(outputs["cna"]).name == "tiny_copykat_CNA_final_results_bin_by_cell.txt"
    assert Path(outputs["manifest"]).exists()

