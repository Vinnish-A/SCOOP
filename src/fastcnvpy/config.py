from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FastCNVConfig:
    """Parameters matching fastCNV's CNVCalling defaults."""

    scale_on_reference_label: bool = True
    threshold_percentile: float = 0.01
    window_size: int = 150
    window_step: int = 10
    top_n_genes: int = 7000
    cluster_k: int | None = None
    cluster_h: float | None = None
    merge_cnv: bool = True
    merge_threshold: float = 0.98
    classification_peaks: tuple[float, float, float] = (-0.1, 0.0, 0.1)

    def validate(self) -> None:
        if not 0 <= self.threshold_percentile < 0.5:
            raise ValueError("threshold_percentile must be in [0, 0.5)")
        if self.window_size < 1:
            raise ValueError("window_size must be >= 1")
        if self.window_step < 1:
            raise ValueError("window_step must be >= 1")
        if self.top_n_genes < 1:
            raise ValueError("top_n_genes must be >= 1")
        if self.cluster_k is not None and self.cluster_k < 1:
            raise ValueError("cluster_k must be >= 1")
        if not 0 <= self.merge_threshold <= 1:
            raise ValueError("merge_threshold must be in [0, 1]")
        if len(self.classification_peaks) != 3:
            raise ValueError("classification_peaks must contain 3 values")

