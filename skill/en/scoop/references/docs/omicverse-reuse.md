# 03. OmicVerse Reuse Strategy

OmicVerse is currently a unified Python multi-omics platform covering bulk, single-cell, spatial, visualization, model-based analysis, and AI-assisted workflows. Its value lies in its many AnnData-native facilities and wrappers, rather than replacing this SOP’s decision logic.

## Recommended reuse

### I/O and data ingestion

Reusable:

- `omicverse.io.read_h5ad`
- `omicverse.io.read_10x_h5`
- `omicverse.io.read_10x_mtx`
- `omicverse.io.read_visium_hd`
- `omicverse.io.read_visium_hd_bin`
- `omicverse.io.read_visium_hd_seg`
- `omicverse.io.read_xenium`
- `omicverse.io.save`

These functions can reduce glue code for input formats from different platforms, but imported data must still undergo this SOP’s schema validation.

### Marker support

Reusable:

- `omicverse.single.find_markers`
- `omicverse.single.get_markers`

`find_markers` supports methods including COSG, wilcoxon, and logreg. COSG can provide fast marker discovery for quickly inspecting large objects; formal evidence tables must still be written externally and enter annotation fusion.

### NMF validation

Reusable:

- `omicverse.single.cNMF`

This SOP uses the FastCNMF profile by default. Invoke the OmicVerse cNMF wrapper only when programme stability is insufficient, a key conclusion depends on programmes, or publication-grade validation is needed.

### CCC validation

Reusable:

- `omicverse.single.run_cellphonedb_v5`
- `omicverse.single.run_liana`
- `omicverse.single.format_liana_results`
- OmicVerse CCC visualization functions

FastCCC is the default initial screen. CellPhoneDB v5 / LIANA is used only for multimers, LRIs central to a mechanism, or ranking robustness validation.

### Report drafts

Reusable:

- `omicverse.single.generate_scRNA_report`

It can generate MultiQC-style report drafts, but cannot replace this SOP’s audit logs and module reports.

## Not recommended as defaults

Not included by default:

- `omicverse.single.lazy`
- CellVote / GPTCelltype as final annotators
- trajectory zoo
- velocity zoo
- drug response
- CNV zoo
- perturbation zoo
- foundation model zoo
- OmicVerse full preprocessing pipeline

These features are not bad tools, but they broaden the SOP’s conceptual scope and complicate routine project delivery. Run them as extended analyses only when a specific biological question requires them.
