from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FastCopyKatConfig:
    """Execution parameters for a CopyKAT-compatible Python run."""

    id_type: str = "S"
    genome: str = "hg20"
    min_gene_per_cell: int = 200
    min_gene_per_chromosome: int = 5
    low_detection_rate: float = 0.05
    upper_detection_rate: float = 0.10
    window_size: int = 25
    bin_size: int = 220_000
    segmentation_threshold: float = 0.10
    min_cluster_cells: int = 5
    max_baseline_clusters: int = 6
    distance: str = "euclidean"
    prediction_score_quantile: float = 0.35
    prediction_mad_multiplier: float = 3.0
    chromosome_rescue_mad_multiplier: float = 3.0
    chromosome_rescue_top_n: int = 6
    random_state: int = 20260615

    def validate(self) -> None:
        if self.min_gene_per_cell < 1:
            raise ValueError("min_gene_per_cell must be >= 1")
        if self.min_gene_per_chromosome < 1:
            raise ValueError("min_gene_per_chromosome must be >= 1")
        if not 0 <= self.low_detection_rate <= 1:
            raise ValueError("low_detection_rate must be in [0, 1]")
        if not 0 <= self.upper_detection_rate <= 1:
            raise ValueError("upper_detection_rate must be in [0, 1]")
        if self.window_size < 3:
            raise ValueError("window_size must be >= 3")
        if self.bin_size < 1:
            raise ValueError("bin_size must be >= 1")
        if self.min_cluster_cells < 1:
            raise ValueError("min_cluster_cells must be >= 1")
        if self.max_baseline_clusters < 2:
            raise ValueError("max_baseline_clusters must be >= 2")
        if self.distance not in {"euclidean", "correlation"}:
            raise ValueError("distance must be 'euclidean' or 'correlation'")
        if self.prediction_mad_multiplier <= 0:
            raise ValueError("prediction_mad_multiplier must be > 0")
        if self.chromosome_rescue_mad_multiplier <= 0:
            raise ValueError("chromosome_rescue_mad_multiplier must be > 0")
        if self.chromosome_rescue_top_n < 1:
            raise ValueError("chromosome_rescue_top_n must be >= 1")
