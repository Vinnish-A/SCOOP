from __future__ import annotations


def require_omicverse():
    """Import OmicVerse from the locked SCOOP runtime."""
    try:
        import omicverse as ov

        return ov
    except Exception as exc:
        raise ImportError("OmicVerse is unavailable in the managed SCOOP runtime; reinstall SCOOP.") from exc


def omicverse_available() -> bool:
    try:
        require_omicverse()
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
