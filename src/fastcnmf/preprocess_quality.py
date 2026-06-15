from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def _as_dense(x) -> np.ndarray:
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x, dtype=np.float64)


def _stream_matrix_stats(reference, candidate, *, chunk_size: int, ref_var_idx=slice(None), cand_var_idx=slice(None)) -> dict[str, float]:
    n_obs = reference.n_obs
    dot = 0.0
    ref_sq = 0.0
    cand_sq = 0.0
    diff_sq = 0.0
    abs_sum = 0.0
    total = 0
    ref_sum = 0.0
    cand_sum = 0.0
    ref_mean_sq_sum = 0.0
    cand_mean_sq_sum = 0.0
    cross_sum = 0.0

    # First pass for means.
    for start in range(0, n_obs, chunk_size):
        end = min(start + chunk_size, n_obs)
        ref = _as_dense(reference.X[start:end, ref_var_idx])
        cand = _as_dense(candidate.X[start:end, cand_var_idx])
        ref_sum += float(ref.sum())
        cand_sum += float(cand.sum())
        total += int(ref.size)
    ref_mean = ref_sum / total
    cand_mean = cand_sum / total

    for start in range(0, n_obs, chunk_size):
        end = min(start + chunk_size, n_obs)
        ref = _as_dense(reference.X[start:end, ref_var_idx])
        cand = _as_dense(candidate.X[start:end, cand_var_idx])
        diff = ref - cand
        dot += float(np.sum(ref * cand))
        ref_sq += float(np.sum(ref * ref))
        cand_sq += float(np.sum(cand * cand))
        diff_sq += float(np.sum(diff * diff))
        abs_sum += float(np.sum(np.abs(diff)))
        ref_centered = ref - ref_mean
        cand_centered = cand - cand_mean
        ref_mean_sq_sum += float(np.sum(ref_centered * ref_centered))
        cand_mean_sq_sum += float(np.sum(cand_centered * cand_centered))
        cross_sum += float(np.sum(ref_centered * cand_centered))

    cosine = dot / max(np.sqrt(ref_sq) * np.sqrt(cand_sq), 1e-12)
    pearson = cross_sum / max(np.sqrt(ref_mean_sq_sum) * np.sqrt(cand_mean_sq_sum), 1e-12)
    return {
        "cosine": float(cosine),
        "pearson": float(pearson),
        "rmse": float(np.sqrt(diff_sq / total)),
        "mean_abs": float(abs_sum / total),
        "n_values": int(total),
    }


def compare_preprocess_outputs(
    *,
    reference_h5ad: Path,
    candidate_h5ad: Path,
    output_json: Path,
    chunk_size: int = 2000,
) -> dict:
    import anndata as ad

    ref = ad.read_h5ad(reference_h5ad, backed="r")
    cand = ad.read_h5ad(candidate_h5ad, backed="r")
    try:
        common_obs = ref.obs_names.intersection(cand.obs_names)
        common_vars = ref.var_names.intersection(cand.var_names)
        if len(common_obs) == 0 or len(common_vars) == 0:
            raise ValueError("reference and candidate have no common observations or variables")

        if not np.array_equal(ref.obs_names.to_numpy(), cand.obs_names.to_numpy()):
            raise ValueError("reference and candidate obs order must match for backed chunked comparison")
        ref_var_idx = [ref.var_names.get_loc(var) for var in common_vars]
        cand_var_idx = [cand.var_names.get_loc(var) for var in common_vars]
        if len(common_vars) == ref.n_vars and np.array_equal(ref.var_names.to_numpy(), common_vars.to_numpy()):
            ref_var_idx = slice(None)
        if len(common_vars) == cand.n_vars and np.array_equal(cand.var_names.to_numpy(), common_vars.to_numpy()):
            cand_var_idx = slice(None)
        matrix = _stream_matrix_stats(ref, cand, chunk_size=chunk_size, ref_var_idx=ref_var_idx, cand_var_idx=cand_var_idx)

        result = {
            "reference_h5ad": str(reference_h5ad),
            "candidate_h5ad": str(candidate_h5ad),
            "reference_shape": [int(ref.n_obs), int(ref.n_vars)],
            "candidate_shape": [int(cand.n_obs), int(cand.n_vars)],
            "common_obs": int(len(common_obs)),
            "common_vars": int(len(common_vars)),
            "var_jaccard": float(len(common_vars) / len(ref.var_names.union(cand.var_names))),
            "obs_jaccard": float(len(common_obs) / len(ref.obs_names.union(cand.obs_names))),
            "matrix": matrix,
            "passes_95pct_input_gate": bool(matrix["cosine"] >= 0.95 and matrix["pearson"] >= 0.95),
        }
    finally:
        ref.file.close()
        cand.file.close()

    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result
