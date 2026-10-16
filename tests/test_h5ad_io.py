
import anndata as ad
import numpy as np
import pytest
from scipy import sparse

from scsp_agent_sop.h5ad_io import H5ADReader
from scsp_agent_sop.qc import inherit_qc


@pytest.mark.parametrize('encoding', ['csr', 'csc', 'dense'])
def test_selected_reads_preserve_counts_order_duplicates_and_empty(tmp_path, encoding):
    matrix = np.array([[0, 2, 0, 4], [3, 0, 5, 0], [0, 0, 0, 0], [7, 8, 9, 0]])
    stored = matrix if encoding == 'dense' else getattr(sparse, f'{encoding}_matrix')(matrix)
    obj = ad.AnnData(stored)
    obj.layers['counts'] = stored.copy()
    path = tmp_path / 'matrix.h5ad'
    obj.write_h5ad(path)
    with H5ADReader(path) as reader:
        for rows, cols in (([3, 0, 3, 2], [2, 0, 2]), ([], [1, 0]), ([2, 1], []), ([True, False, True, False], slice(None))):
            actual = reader.read_rows(rows, cols, layer='counts')
            actual = actual.toarray() if sparse.issparse(actual) else actual
            expected = matrix[np.asarray(rows, dtype=bool if rows and isinstance(rows[0], bool) else int)][:, cols]
            np.testing.assert_array_equal(actual, expected)
        np.testing.assert_array_equal(reader.read_rows(slice(1, 4, 2)).toarray() if encoding != 'dense' else reader.read_rows(slice(1, 4, 2)), matrix[1:4:2])
        with pytest.raises(IndexError):
            reader.read_rows([4])


def test_inherited_qc_verifies_exact_counts_and_retained_membership(tmp_path):
    parent = ad.AnnData(sparse.csr_matrix(np.arange(24).reshape(6, 4)))
    parent.layers['counts'] = parent.X.copy()
    for key, value in {'sample_id': 's1', 'qc_pass': True, 'qc_class': 'pass', 'final_use': True,
                       'doublet_score': .1, 'doublet_call_scrublet': False, 'doublet_call': 'singlet'}.items():
        parent.obs[key] = value
    parent.obs.loc['0', 'final_use'] = False
    path = tmp_path / 'parent.h5ad'
    parent.write_h5ad(path)
    child = parent[[4, 2]].copy()
    child.obs['doublet_score'] = 99.
    record = inherit_qc(child, path, chunk_size=1)
    assert record['mode'] == 'inherited'
    np.testing.assert_array_equal(child.obs.doublet_score, [.1, .1])
    child.layers['counts'][0, 1] += 1
    with pytest.raises(ValueError, match='exact parent counts'):
        inherit_qc(child, path)
    with pytest.raises(ValueError, match='final_use'):
        inherit_qc(parent[[0]].copy(), path)


def test_managed_dependency_contract_supports_native_backed_indexing(tmp_path):
    obj = ad.AnnData(sparse.csr_matrix([[0, 2], [3, 0], [1, 4]]))
    path = tmp_path / 'backed.h5ad'
    obj.write_h5ad(path)
    backed = ad.read_h5ad(path, backed='r')
    try:
        np.testing.assert_array_equal(backed[[2, 0, 2], :].X.toarray(), [[1, 4], [0, 2], [1, 4]])
    finally:
        backed.file.close()
