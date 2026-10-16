import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from scsp_agent_sop.qc import compute_basic_qc, run_scrublet_per_sample


def test_real_scrublet_parallel_matches_serial_and_resumes(tmp_path):
    rng = np.random.default_rng(32)
    counts = sparse.csr_matrix(rng.poisson(rng.gamma(1., 2., size=(600, 500))))
    obj = ad.AnnData(counts.copy(), obs=pd.DataFrame({'sample_id': ['a'] * 300 + ['b'] * 300}, index=[f'c{i}' for i in range(600)]))
    obj.layers['counts'] = counts
    obj.obs['final_use'] = True
    compute_basic_qc(obj)
    serial, parallel = obj.copy(), obj.copy()
    run_scrublet_per_sample(serial, workers=1)
    first = run_scrublet_per_sample(parallel, workers=2, checkpoint_dir=tmp_path)
    np.testing.assert_array_equal(serial.obs.doublet_score, parallel.obs.doublet_score)
    np.testing.assert_array_equal(serial.obs.doublet_call, parallel.obs.doublet_call)
    resumed = run_scrublet_per_sample(obj, workers=2, checkpoint_dir=tmp_path)
    assert resumed.resumed.all()
    assert len(list(tmp_path.glob('*.npz'))) == 2
    np.testing.assert_array_equal(obj.obs.doublet_score, parallel.obs.doublet_score)
    obj.layers['counts'][0, 0] += 1
    changed = run_scrublet_per_sample(obj, workers=2, checkpoint_dir=tmp_path)
    assert changed.resumed.sum() == 1
    assert len(list(tmp_path.glob('*.npz'))) == 3
    assert set(first.backend) == {'scanpy_scrublet'}
