# Annotation evidence


Inputs:

- `artifacts/adata_programs.h5ad`.
- Final clusters and NMF programmes.

Script:

```bash
python scripts/04_annotation_markers.py --config runs/<run_id>/config/run.yaml
```

OmicVerse fast marker option:

```bash
python scripts/04_annotation_markers.py --config runs/<run_id>/config/run.yaml --use-omicverse-cosg
```

The OmicVerse wrapper is used only to export marker evidence. Final labels still require evidence fusion.

H5AD writes:

- No final labels unless supplied by user/Agent after review.

Sidecar outputs:

- `tables/cluster_markers.parquet`.
- `tables/annotation_evidence_template.tsv`.

Review triggers:

- Reference and marker evidence disagree.
- New cell type claim.
- Low-confidence annotation used in downstream conclusions.
