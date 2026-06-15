from __future__ import annotations

from dataclasses import dataclass


DTYPE_BYTES = {
    "float16": 2,
    "bfloat16": 2,
    "float32": 4,
    "float64": 8,
}


@dataclass(frozen=True)
class NMFMemoryEstimate:
    """Approximate memory footprint for one NMF factorization task."""

    observations: int
    genes: int
    components: int
    dtype: str
    replicate_batch_size: int
    estimated_bytes: int

    @property
    def estimated_gib(self) -> float:
        return self.estimated_bytes / 1024**3


def estimate_nmf_memory(
    observations: int,
    genes: int,
    components: int,
    dtype: str = "float32",
    replicate_batch_size: int = 1,
    temporary_multiplier: float = 3.0,
) -> NMFMemoryEstimate:
    """Estimate memory required by a dense NMF task.

    The estimate intentionally errs high: it includes X, W, H, and temporary
    buffers. GPU backends use this before choosing chunk and replicate batch
    sizes.
    """

    if dtype not in DTYPE_BYTES:
        raise ValueError(f"unsupported dtype: {dtype}")
    if min(observations, genes, components, replicate_batch_size) < 1:
        raise ValueError("matrix dimensions and replicate_batch_size must be positive")

    itemsize = DTYPE_BYTES[dtype]
    x = observations * genes * itemsize
    w = observations * components * itemsize * replicate_batch_size
    h = components * genes * itemsize * replicate_batch_size
    estimated = int((x + w + h) * temporary_multiplier)
    return NMFMemoryEstimate(
        observations=observations,
        genes=genes,
        components=components,
        dtype=dtype,
        replicate_batch_size=replicate_batch_size,
        estimated_bytes=estimated,
    )


def choose_replicate_batch_size(
    observations: int,
    genes: int,
    components: int,
    available_bytes: int,
    dtype: str = "float32",
    safety_fraction: float = 0.8,
    max_batch_size: int = 8,
) -> int:
    """Choose the largest replicate batch that fits the memory budget."""

    if available_bytes <= 0:
        raise ValueError("available_bytes must be positive")
    if not 0 < safety_fraction <= 1:
        raise ValueError("safety_fraction must be in (0, 1]")

    budget = int(available_bytes * safety_fraction)
    for batch_size in range(max_batch_size, 0, -1):
        estimate = estimate_nmf_memory(
            observations=observations,
            genes=genes,
            components=components,
            dtype=dtype,
            replicate_batch_size=batch_size,
        )
        if estimate.estimated_bytes <= budget:
            return batch_size
    return 1

