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
    # Upper bound, in stored entries, on one coalesced read.
    max_span = 1 << 26

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
        pieces, indices = [], []
        indptr = np.zeros(len(unique) + 1, dtype=np.int64)
        if len(unique):
            ptr = node["indptr"][int(unique[0]):int(unique[-1]) + 2].astype(np.int64)
            starts, ends = ptr[unique - unique[0]], ptr[unique - unique[0] + 1]
            indptr[1:] = np.cumsum(ends - starts)
            # Every hyperslab re-reads its boundary chunks, so isolated rows cost
            # a chunk each. Read through gaps shorter than a chunk instead, and
            # cap the span so a dense selection never loads the whole matrix.
            chunks = node["data"].chunks
            max_gap = chunks[0] if chunks else 1 << 16
            max_span = max(max_gap, self.max_span)
            breaks = np.r_[np.flatnonzero(starts[1:] - ends[:-1] > max_gap) + 1, len(unique)]
            head = 0
            while head < len(unique):
                stop = int(breaks[np.searchsorted(breaks, head, "right")])
                stop = max(head + 1, min(stop, int(np.searchsorted(ends, starts[head] + max_span, "right"))))
                a, b = int(starts[head]), int(ends[stop - 1])
                data, idx = node["data"][a:b], node["indices"][a:b]
                if indptr[stop] - indptr[head] != b - a:
                    lengths = ends[head:stop] - starts[head:stop]
                    take = np.repeat(starts[head:stop] - a - (indptr[head:stop] - indptr[head]), lengths)
                    take += np.arange(indptr[stop] - indptr[head])
                    data, idx = data[take], idx[take]
                pieces.append(data)
                indices.append(idx)
                head = stop
        data = np.concatenate(pieces) if pieces else np.empty(0, dtype=node["data"].dtype)
        idx = np.concatenate(indices) if indices else np.empty(0, dtype=node["indices"].dtype)
        minor_size = self.n_vars if encoding == "csr_matrix" else self.n_obs
        matrix = sparse.csr_matrix((data, idx, indptr), shape=(len(unique), minor_size))
        matrix = matrix[inverse][:, minor]
        return matrix if encoding == "csr_matrix" else matrix.T.tocsr()
