# References Checked for This Template

The project templates and wrappers were designed around the following current documentation and method facts:

- OmicVerse v2 is documented as a unified Python framework covering bulk RNA-seq, single-cell, spatial transcriptomics, visualisation, model-based analysis and AI-assisted workflows.
- OmicVerse API lists reusable I/O functions such as `read_h5ad`, `read_10x_h5`, `read_10x_mtx`, `read_visium_hd`, `read_visium_hd_bin`, `read_visium_hd_seg`, `read_xenium` and `save`.
- OmicVerse API lists preprocessing utilities including `anndata_to_GPU`, `anndata_to_CPU`, `scrublet`, PCA, neighbours, UMAP and Leiden. In this SOP, these are used as infrastructure or fallback wrappers, not as an all-in-one pipeline.
- `omicverse.single.run_cellphonedb_v5` accepts an AnnData, a CellPhoneDB zip file path, cell type key and parameters such as `iterations`, `threshold`, `pvalue` and `threads`. The wrapper writes raw CellPhoneDB result objects and a communication AnnData to `adata.uns`, so this project immediately exports those results to external files and removes heavy objects from the H5AD.
- `omicverse.single.run_liana` accepts an AnnData, `groupby`, `method='rank_aggregate'`, `key_added`, `inplace` and additional kwargs. This is used as a validation wrapper rather than as the primary CCC screen.
- `omicverse.single.cNMF` is a consensus NMF workflow wrapper with candidate `components`, `n_iter`, `output_dir`, `use_gpu` and `gpu_id`; this project uses it as a validation path for unstable programmes, not as the default NMF step.
- FastCCC is a scalable, permutation-free framework for CCC detection. In this SOP it is the primary CCC screening tool, while multimeric or mechanism-critical interactions require CellPhoneDB/LIANA validation.
