from __future__ import annotations

import sys
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _submodule_path() -> Path:
    return _repo_root() / "third_party" / "omicverse"


def require_omicverse(*, source: str = "auto", allow_submodule: bool = True):
    """Lazy import OmicVerse from an installation or the pinned submodule."""
    if source not in {"auto", "installed", "submodule"}:
        raise ValueError("source must be one of: auto, installed, submodule")
    if source in {"auto", "installed"}:
        try:
            import omicverse as ov

            return ov
        except Exception as installed_exc:
            if source == "installed" or not allow_submodule:
                raise ImportError("OmicVerse is not installed.") from installed_exc
    path = _submodule_path()
    if source in {"auto", "submodule"} and allow_submodule and path.exists():
        path_str = str(path)
        if path_str not in sys.path:
            sys.path.insert(0, path_str)
        try:
            import omicverse as ov

            return ov
        except Exception as submodule_exc:
            raise ImportError(f"OmicVerse submodule is present but could not be imported from {path}.") from submodule_exc
    raise ImportError("OmicVerse is unavailable. Install it or initialize third_party/omicverse.")


def omicverse_available(*, allow_submodule: bool = True) -> bool:
    try:
        require_omicverse(allow_submodule=allow_submodule)
        return True
    except Exception:
        return False


def read_h5ad(path: str):
    ov = require_omicverse()
    return ov.io.read_h5ad(path)


def save(adata, path: str):
    ov = require_omicverse()
    return ov.io.save(adata, path)


def generate_report(adata, output_path: str, species: str = "human", sample_key: str | None = None):
    ov = require_omicverse()
    return ov.single.generate_scRNA_report(adata, output_path=output_path, species=species, sample_key=sample_key)


def run_cnmf_validation(adata, components, output_dir: str, n_iter: int = 100, use_gpu: bool = True):
    ov = require_omicverse()
    return ov.single.cNMF(
        adata,
        components=components,
        n_iter=n_iter,
        output_dir=output_dir,
        use_gpu=use_gpu,
        name="omicverse_cnmf_validation",
    )


__all__ = [
    "generate_report",
    "omicverse_available",
    "read_h5ad",
    "require_omicverse",
    "run_cnmf_validation",
    "save",
]
