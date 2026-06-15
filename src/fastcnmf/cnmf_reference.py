from __future__ import annotations

import json
from pathlib import Path


def run_cnmf_reference_preprocess(
    *,
    input_h5ad: Path,
    output_prefix: Path,
    sample_key: str = "sample_id",
    n_top_genes: int = 3000,
    librarysize_targetsum: float = 1e4,
    theta: float = 1.0,
    max_iter_harmony: int = 20,
    seed: int = 20260614,
) -> dict:
    """Run the upstream cNMF preprocessing contract from a cold H5AD input."""

    import anndata as ad
    from cnmf.preprocess import Preprocess

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(input_h5ad)
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    if sample_key not in adata.obs.columns:
        raise ValueError(f"sample_key {sample_key!r} is not present in obs")

    pre = Preprocess(random_seed=seed)
    corrected, tp10k, hvgs = pre.preprocess_for_cnmf(
        adata,
        harmony_vars=sample_key,
        n_top_rna_genes=n_top_genes,
        librarysize_targetsum=librarysize_targetsum,
        makeplots=False,
        theta=theta,
        save_output_base=str(output_prefix),
        max_iter_harmony=max_iter_harmony,
    )
    corrected_path = Path(str(output_prefix) + ".Corrected.HVG.Varnorm.h5ad")
    tp10k_path = Path(str(output_prefix) + ".TP10K.h5ad")
    hvg_path = Path(str(output_prefix) + ".Corrected.HVGs.txt")
    manifest = {
        "engine": "cnmf.preprocess.Preprocess.preprocess_for_cnmf",
        "input_h5ad": str(input_h5ad),
        "output_prefix": str(output_prefix),
        "corrected_h5ad": str(corrected_path),
        "tp10k_h5ad": str(tp10k_path),
        "hvg_txt": str(hvg_path),
        "sample_key": sample_key,
        "n_obs": int(corrected.n_obs),
        "n_vars": int(corrected.n_vars),
        "n_hvgs": int(len(hvgs)),
        "n_top_genes": int(n_top_genes),
        "librarysize_targetsum": float(librarysize_targetsum),
        "theta": float(theta),
        "max_iter_harmony": int(max_iter_harmony),
        "seed": int(seed),
    }
    manifest_path = output_prefix.parent / "cnmf_reference_preprocess_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
