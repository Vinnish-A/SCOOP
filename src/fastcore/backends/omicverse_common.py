from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping


def require_omicverse():
    try:
        import omicverse as ov
    except Exception as exc:  # pragma: no cover - optional dependency.
        raise RuntimeError("OmicVerse is not available for this FastCore backend.") from exc
    return ov


def unsupported_backend_result(backend: str, cfg: Mapping[str, Any], run_root: str | Path) -> dict[str, Any]:
    raise RuntimeError(
        f"{backend} was selected, but its executable OmicVerse adapter is not enabled in this PR. "
        "Use scanpy_legacy fallback or install and validate the OmicVerse adapter before selecting this backend."
    )


def map_standard_core_keys(adata) -> None:
    """Map common OmicVerse keys back to SCOOP's stable downstream schema."""
    if "X_pca" in adata.obsm and "X_pca_biology" not in adata.obsm:
        adata.obsm["X_pca_biology"] = adata.obsm["X_pca"].copy()
    if "X_pca" in adata.obsm and "X_pca_identity_prebatch" not in adata.obsm:
        adata.obsm["X_pca_identity_prebatch"] = adata.obsm["X_pca"].copy()
    if "X_umap" in adata.obsm:
        adata.obsm.setdefault("X_umap_biology", adata.obsm["X_umap"].copy())
        adata.obsm.setdefault("X_umap_identity", adata.obsm["X_umap"].copy())
    if "connectivities" in adata.obsp:
        adata.obsp.setdefault("connectivities_identity", adata.obsp["connectivities"].copy())
    if "distances" in adata.obsp:
        adata.obsp.setdefault("distances_identity", adata.obsp["distances"].copy())
