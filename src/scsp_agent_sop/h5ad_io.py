"""Read H5AD axes and selected matrix rows without eagerly loading layers.

Uses the on-disk CSR/CSC contract rather than AnnData backed sparse indexing.
Files are read-only; row order, repeated requests and sparse duplicates survive.
"""
from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
from anndata.io import read_elem
from scipy import sparse


def _positions(index, size):
    if isinstance(index, slice):
        return np.arange(*index.indices(size), dtype=np.intp)
    values = np.asarray(index)
    if values.dtype == bool:
        if values.shape != (size,):
            raise IndexError("boolean index must match the matrix axis")
        return np.flatnonzero(values)
    if values.size == 0:
        return np.empty(0, dtype=np.intp)
    if values.ndim != 1 or values.dtype.kind not in "iu":
        raise IndexError("matrix indices must be a one-dimensional integer or boolean array")
    values = values.astype(np.intp, copy=True)
    values[values < 0] += size
    if np.any((values < 0) | (values >= size)):
        raise IndexError("matrix index out of bounds")
    return values


class H5ADReader:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.file = h5py.File(self.path, "r")
        try:
            self.obs = read_elem(self.file["obs"])
            self.var = read_elem(self.file["var"])
        except Exception:
            self.file.close()
            raise
        self.obs_names, self.var_names = self.obs.index, self.var.index
        self.n_obs, self.n_vars = len(self.obs), len(self.var)
        self.shape = (self.n_obs, self.n_vars)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        self.file.close()

    def read_rows(self, rows=slice(None), columns=slice(None), *, layer=None):
        node = self.file["X" if layer is None else f"layers/{layer}"]
        row = _positions(rows, self.n_obs)
        col = _positions(columns, self.n_vars)
        if isinstance(node, h5py.Dataset):
            # h5py fancy indices must be unique and sorted.
            unique, inverse = np.unique(row, return_inverse=True)
            return node[unique, :][inverse][:, col]
        encoding = node.attrs["encoding-type"]
        if isinstance(encoding, bytes):
            encoding = encoding.decode()
        if encoding not in {"csr_matrix", "csc_matrix"}:
            raise ValueError(f"unsupported H5AD matrix encoding: {encoding}")
        major = row if encoding == "csr_matrix" else col
        minor = col if encoding == "csr_matrix" else row
        unique, inverse = np.unique(major, return_inverse=True)
        pieces, indices, indptr = [], [], [0]
        # Read contiguous requests together, without reading unrelated layers.
        for run in np.split(unique, np.flatnonzero(np.diff(unique) != 1) + 1):
            if not len(run):
                continue
            ptr = node["indptr"][int(run[0]):int(run[-1]) + 2]
            start, end = int(ptr[0]), int(ptr[-1])
            pieces.append(node["data"][start:end])
            indices.append(node["indices"][start:end])
            indptr.extend((ptr[1:] - start + indptr[-1]).tolist())
        data = np.concatenate(pieces) if pieces else np.empty(0, dtype=node["data"].dtype)
        idx = np.concatenate(indices) if indices else np.empty(0, dtype=np.int32)
        minor_size = self.n_vars if encoding == "csr_matrix" else self.n_obs
        matrix = sparse.csr_matrix((data, idx, indptr), shape=(len(unique), minor_size))
        matrix = matrix[inverse][:, minor]
        return matrix if encoding == "csr_matrix" else matrix.T.tocsr()
