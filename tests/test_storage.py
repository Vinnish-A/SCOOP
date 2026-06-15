from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd

from scsp_agent_sop.storage import prune_h5ad


def test_prune_h5ad_minimal_schema():
    a = ad.AnnData(np.ones((2, 2)))
    a.obs["sample_id"] = ["s1", "s1"]
    a.obs["tmp"] = [1, 2]
    a.var["gene_symbol"] = ["A", "B"]
    a.var["tmp"] = [1, 2]
    a.uns["file_registry"] = {}
    a.uns["big_tmp"] = {"x": list(range(10))}
    schema = {"obs_keep": {"identity": ["sample_id"]}, "var_keep": ["gene_symbol"], "layers_keep": [], "obsm_keep": [], "obsp_keep": [], "uns_keep": ["file_registry"]}
    prune_h5ad(a, schema)
    assert "sample_id" in a.obs
    assert "tmp" not in a.obs
    assert "gene_symbol" in a.var
    assert "big_tmp" not in a.uns
