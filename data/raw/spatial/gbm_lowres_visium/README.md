# Low-resolution Visium raw fixtures

This directory is reserved for de-identified low-resolution Visium samples used
by raw spatial transcriptomics input tests. Real sample folders are intentionally
not tracked in git.

Each sample keeps the standard 10x matrix inputs and only the low-resolution
spatial image files needed by common Visium readers:

- `filtered_feature_bc_matrix.h5`
- `filtered_feature_bc_matrix/{barcodes.tsv.gz,features.tsv.gz,matrix.mtx.gz}`
- `spatial/tissue_lowres_image.png`
- `spatial/scalefactors_json.json`
- `spatial/tissue_positions.csv`

