from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class HarmonyConfig:
    """Harmony preprocessing contract for cNMF-compatible mode."""

    batch_key: str = "sample_id"
    n_top_genes: int = 3000
    librarysize_targetsum: float = 1e4
    theta: float = 1.0
    max_iter_harmony: int = 20
    harmonypy_version: str = "0.2.0"


@dataclass(frozen=True)
class NMFConfig:
    """NMF search grid and execution parameters."""

    k_values: tuple[int, ...] = (6, 8, 10, 12)
    n_iter: int = 20
    max_nmf_iter: int = 50
    seed: int = 20260614
    dtype: str = "float64"
    mode: str = "exact-compatible"


@dataclass(frozen=True)
class ResourceConfig:
    """Runtime resource controls."""

    workers: int = 4
    blas_threads_per_worker: int = 1
    gpu_enabled: bool = False
    gpu_memory_safety_fraction: float = 0.8
    fallback_to_cpu: bool = True


@dataclass(frozen=True)
class FastCNMFConfig:
    """Top-level FastCNMF run configuration."""

    input_h5ad: Path
    output_dir: Path
    run_name: str = "fastcnmf"
    harmony: HarmonyConfig = field(default_factory=HarmonyConfig)
    nmf: NMFConfig = field(default_factory=NMFConfig)
    resources: ResourceConfig = field(default_factory=ResourceConfig)

    def validate(self) -> None:
        if not self.input_h5ad:
            raise ValueError("input_h5ad is required")
        if not self.output_dir:
            raise ValueError("output_dir is required")
        if not self.nmf.k_values:
            raise ValueError("at least one k value is required")
        if self.nmf.n_iter < 1:
            raise ValueError("n_iter must be >= 1")
        if self.resources.workers < 1:
            raise ValueError("workers must be >= 1")
        if not 0 < self.resources.gpu_memory_safety_fraction <= 1:
            raise ValueError("gpu_memory_safety_fraction must be in (0, 1]")
