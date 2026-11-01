"""Per-sample feature measurement (``scoop.feature-measurement.v1``).

A gene that a sample's capture did not quantify is unmeasured there, not a
zero. ``FeatureMeasurement`` records ``measured[sample, gene]`` for one var
axis so consumers can restrict per-gene statistics to measuring cells. Cell
masks are derived through sample codes; the object never stores a
cell x gene matrix.

On disk it is an ``.npz`` with the boolean matrix bit-packed along genes.
An AnnData declares it at ``uns["file_registry"]["artifacts"]["feature_measurement"]``
(``storage.register_file(..., category="artifacts", schema=FORMAT)``).
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

FORMAT = "scoop.feature-measurement.v1"
REGISTRY_KEY = "feature_measurement"


def var_sha256(var_names) -> str:
    return hashlib.sha256("\n".join(map(str, var_names)).encode("utf-8")).hexdigest()


@dataclass(frozen=True, eq=False)
class FeatureMeasurement:
    sample_key: str
    samples: pd.Index
    var_sha256: str
    measured: np.ndarray
    reasons: np.ndarray | None = None
    reason_legend: Mapping[int, str] | None = None

    def __post_init__(self):
        samples = pd.Index(np.asarray(self.samples).astype(str))
        if not samples.is_unique:
            raise ValueError("Feature measurement samples must be unique")
        measured = np.asarray(self.measured, dtype=bool)
        if measured.ndim != 2 or measured.shape[0] != len(samples):
            raise ValueError(f"measured must be samples x genes; got {measured.shape} for {len(samples)} samples")
        if self.reasons is not None and np.shape(self.reasons) != measured.shape:
            raise ValueError("reasons must have the shape of measured")
        object.__setattr__(self, "samples", samples)
        object.__setattr__(self, "measured", measured)

    @classmethod
    def from_var(cls, sample_key: str, samples, var_names, measured, reasons=None, reason_legend=None) -> "FeatureMeasurement":
        measured = np.asarray(measured, dtype=bool)
        if measured.shape[1:] != (len(var_names),):
            raise ValueError(f"measured has {measured.shape[1:]} genes; var axis has {len(var_names)}")
        return cls(sample_key, samples, var_sha256(var_names), measured, reasons, reason_legend)

    @property
    def n_features(self) -> int:
        return self.measured.shape[1]

    def save(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        arrays = {"format": np.array(FORMAT), "sample_key": np.array(self.sample_key),
                  "samples": self.samples.to_numpy().astype(str), "var_sha256": np.array(self.var_sha256),
                  "n_features": np.array(self.n_features), "measured_packbits": np.packbits(self.measured, axis=1)}
        if self.reasons is not None:
            arrays["reasons"] = np.asarray(self.reasons, dtype=np.uint8)
            arrays["reason_legend"] = np.array(json.dumps({str(k): v for k, v in (self.reason_legend or {}).items()}))
        with path.open("wb") as fh:  # a file handle keeps numpy from appending ".npz"
            np.savez_compressed(fh, **arrays)
        return path

    @classmethod
    def load(cls, path: str | Path) -> "FeatureMeasurement":
        with np.load(path, allow_pickle=False) as z:
            if "format" not in z.files or str(z["format"]) != FORMAT:
                raise ValueError(f"{path} is not a {FORMAT} file")
            measured = np.unpackbits(z["measured_packbits"], axis=1, count=int(z["n_features"])).astype(bool)
            reasons = z["reasons"] if "reasons" in z.files else None
            legend = {int(k): v for k, v in json.loads(str(z["reason_legend"])).items()} if "reason_legend" in z.files else None
            return cls(str(z["sample_key"]), z["samples"], str(z["var_sha256"]), measured, reasons, legend)

    def sample_codes(self, obs_samples) -> np.ndarray:
        values = pd.Index(pd.Series(obs_samples).astype(str))
        codes = self.samples.get_indexer(values)
        if (codes < 0).any():
            missing = sorted(set(values[codes < 0]))
            raise ValueError(f"Feature measurement does not cover {int((codes < 0).sum())} cells "
                             f"from samples {missing[:5]}{' ...' if len(missing) > 5 else ''}")
        return codes

    def measured_in_all(self, obs_samples) -> np.ndarray:
        """Genes measured in every sample present among ``obs_samples``."""
        return self.measured[np.unique(self.sample_codes(obs_samples))].all(axis=0)

    def measured_cells(self, obs_samples, genes=slice(None)) -> np.ndarray:
        """Number of cells that measured each gene."""
        cells_per_sample = np.bincount(self.sample_codes(obs_samples), minlength=len(self.samples))
        return cells_per_sample @ self.measured[:, genes]

    def cell_mask(self, obs_samples, genes) -> np.ndarray:
        """Boolean cells x ``genes`` mask; keep ``genes`` small."""
        return self.measured[:, genes][self.sample_codes(obs_samples)]

    def check(self, adata) -> None:
        if self.sample_key not in adata.obs:
            raise KeyError(f"Feature measurement sample_key {self.sample_key!r} is not in obs")
        if var_sha256(adata.var_names) != self.var_sha256:
            raise ValueError("Feature measurement was built for a different var axis (var_sha256 mismatch)")
        self.sample_codes(adata.obs[self.sample_key])


def feature_measurement_path(adata, path: str | Path | None = None) -> Path | None:
    """Explicit path, else the registered artifact, else None (measurement not declared)."""
    if path is not None:
        return Path(path)
    record = adata.uns.get("file_registry", {}).get("artifacts", {}).get(REGISTRY_KEY)
    return Path(record["path"]) if record else None


def measurement_record(adata, *, honoured: bool) -> dict:
    """Stage receipt entry: ``applied``, ``not_honoured`` (declared, but the method cannot use it) or ``not_declared``."""
    path = feature_measurement_path(adata)
    if path is None:
        return {"status": "not_declared"}
    return {"status": "applied" if honoured else "not_honoured", "path": str(path)}
