from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class DatasetSpec:
    """Dataset entry used by reproducible FastCNMF benchmarks."""

    dataset_id: str
    tier: str
    modality: str
    path: str
    source_kind: str
    n_obs: int
    n_vars: int | None = None
    sample_key: str = "sample_id"
    n_samples: int | None = None
    min_cells_per_sample: int | None = None
    median_cells_per_sample: float | None = None
    max_cells_per_sample: int | None = None
    role: str = "candidate_and_reference"
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class BenchmarkLane:
    """One execution lane in a fair benchmark comparison."""

    lane_id: str
    engine: str
    cold_start: bool
    allow_existing_intermediates: bool
    allowed_hardware_acceleration: tuple[str, ...]
    description: str


@dataclass(frozen=True)
class FairnessPolicy:
    """Benchmark rules that prevent accidental warm-cache comparisons."""

    comparison_name: str = "same_hardware_cold_start"
    reuse_policy: str = "no_cross_lane_artifact_reuse"
    reference_lane: str = "cnmf_optimized"
    candidate_lane: str = "fastcnmf_independent"
    target_speedup: float = 3.0
    min_spectra_cosine: float = 0.95
    min_usage_pearson: float = 0.95
    require_same_input_cells_and_genes: bool = True
    require_same_harmony_batch_key: bool = True


@dataclass(frozen=True)
class FastCNMFBenchmarkManifest:
    """Top-level benchmark manifest for S1/S2 FastCNMF development."""

    schema_version: str
    run_id: str
    generated_from: str
    datasets: tuple[DatasetSpec, ...]
    lanes: tuple[BenchmarkLane, ...]
    fairness: FairnessPolicy = field(default_factory=FairnessPolicy)
    parameters: dict[str, Any] = field(default_factory=dict)
    output_root: str = "tmp/fastcnmf_large_benchmark"
    missing_tiers: tuple[str, ...] = ()

    def to_json(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2), encoding="utf-8")

    @classmethod
    def from_json(cls, path: Path) -> "FastCNMFBenchmarkManifest":
        payload = json.loads(path.read_text(encoding="utf-8"))
        return cls(
            schema_version=payload["schema_version"],
            run_id=payload["run_id"],
            generated_from=payload["generated_from"],
            datasets=tuple(DatasetSpec(**item) for item in payload["datasets"]),
            lanes=tuple(BenchmarkLane(**item) for item in payload["lanes"]),
            fairness=FairnessPolicy(**payload["fairness"]),
            parameters=payload.get("parameters", {}),
            output_root=payload.get("output_root", "tmp/fastcnmf_large_benchmark"),
            missing_tiers=tuple(payload.get("missing_tiers", ())),
        )


def _sample_stats(obs, sample_key: str) -> dict[str, int | float]:
    values = obs[sample_key].astype(str).value_counts()
    return {
        "n_samples": int(values.shape[0]),
        "min_cells_per_sample": int(values.min()),
        "median_cells_per_sample": float(values.median()),
        "max_cells_per_sample": int(values.max()),
    }


def inspect_h5ad_dataset(path: Path, *, dataset_id: str, tier: str, modality: str) -> DatasetSpec:
    import anndata as ad

    adata = ad.read_h5ad(path, backed="r")
    try:
        sample_key = "sample_id" if "sample_id" in adata.obs.columns else "sample"
        stats = _sample_stats(adata.obs, sample_key) if sample_key in adata.obs.columns else {}
        return DatasetSpec(
            dataset_id=dataset_id,
            tier=tier,
            modality=modality,
            path=str(path),
            source_kind="h5ad",
            n_obs=int(adata.n_obs),
            n_vars=int(adata.n_vars),
            sample_key=sample_key,
            **stats,
        )
    finally:
        adata.file.close()


def inspect_visium_root(path: Path, *, dataset_id: str, tier: str) -> DatasetSpec:
    samples = sorted(p for p in path.glob("GBM_*") if p.is_dir())
    counts = []
    for sample in samples:
        h5 = sample / "filtered_feature_bc_matrix.h5"
        if not h5.exists():
            continue
        import h5py

        with h5py.File(h5, "r") as handle:
            shape = handle["matrix/shape"][:]
            counts.append(int(shape[1]))
    return DatasetSpec(
        dataset_id=dataset_id,
        tier=tier,
        modality="spatial_visium",
        path=str(path),
        source_kind="visium_directory",
        n_obs=int(sum(counts)),
        n_vars=None,
        sample_key="sample_id",
        n_samples=len(counts),
        min_cells_per_sample=min(counts) if counts else None,
        median_cells_per_sample=float(sorted(counts)[len(counts) // 2]) if counts else None,
        max_cells_per_sample=max(counts) if counts else None,
        notes=("S1 is incomplete until more GBM Visium samples are available.",),
    )


def default_lanes() -> tuple[BenchmarkLane, ...]:
    return (
        BenchmarkLane(
            lane_id="cnmf_optimized",
            engine="cnmf==1.7.1",
            cold_start=True,
            allow_existing_intermediates=False,
            allowed_hardware_acceleration=("multi_process_cpu", "blas_thread_limits"),
            description="Fair same-hardware cNMF baseline; uses all allowed CPU parallelism from raw input.",
        ),
        BenchmarkLane(
            lane_id="fastcnmf_independent",
            engine="fastcnmf",
            cold_start=True,
            allow_existing_intermediates=False,
            allowed_hardware_acceleration=("multi_process_cpu", "blas_thread_limits", "gpu_if_enabled"),
            description="Independent FastCNMF runner; no reuse of cNMF reference intermediates.",
        ),
    )


def build_default_manifest(root: Path, output_root: Path) -> FastCNMFBenchmarkManifest:
    quick = root / "h5ad/canonical/quick_test"
    spatial = root / "data/raw/spatial/gbm_lowres_visium"
    datasets: list[DatasetSpec] = []

    if spatial.exists():
        datasets.append(inspect_visium_root(spatial, dataset_id="s1_gbm_lowres_visium_available", tier="S1"))

    h5ads = {
        "s2_public_24x3000": quick / "public_O_GSE154795_24samples_3000cells_balanced.h5ad",
        "s2_internal_24x3000": quick / "internal_overall_sim_24samples_3000cells_balanced.h5ad",
    }
    for dataset_id, path in h5ads.items():
        if path.exists():
            datasets.append(inspect_h5ad_dataset(path, dataset_id=dataset_id, tier="S2", modality="single_cell"))

    missing = []
    if not any(d.tier == "S1" and d.n_samples and d.n_samples > 3 for d in datasets):
        missing.append("S1_all_gbm_visium")
    if not any(d.tier == "S2" for d in datasets):
        missing.append("S2_72k_h5ad")

    return FastCNMFBenchmarkManifest(
        schema_version="fastcnmf.benchmark_manifest.v1",
        run_id="fastcnmf_s1_s2_cold_start",
        generated_from=str(root),
        datasets=tuple(datasets),
        lanes=default_lanes(),
        fairness=FairnessPolicy(),
        parameters={
            "harmony_batch_key": "sample_id",
            "n_top_genes": 3000,
            "k_values": [6, 8, 10, 12],
            "n_iter_smoke": 8,
            "n_iter_production": 20,
            "max_nmf_iter_smoke": 50,
            "max_nmf_iter_production": 50,
        },
        output_root=str(output_root),
        missing_tiers=tuple(missing),
    )
