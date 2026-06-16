# OmicVerse Transfer Module

`omicverse_transfer` isolates all SCOOP code that is either derived from
OmicVerse or calls optional external OmicVerse APIs.

Active contents:

- `omicverse_transfer.vendor.omicverse_gpl`: vendored GPL CPU preprocessing,
  HVG, PCA, graph, UMAP, and Leiden subset used by `omicverse_cpu`.
- `omicverse_transfer.core_common`: shared FastCore adapter options, key
  mapping, Harmony bridge, manifest/table writing, and external OmicVerse step
  wrappers.
- `omicverse_transfer.core_cpu`: self-contained vendored CPU backend entry
  point.
- `omicverse_transfer.core_mixed`: optional external OmicVerse CPU/GPU mixed
  backend.
- `omicverse_transfer.core_rust_oom`: optional external OmicVerse Rust/OOM
  backend.
- `omicverse_transfer.external`: central optional `import omicverse` gateway for
  report, marker, cNMF validation, CellPhoneDB, and LIANA wrappers.

Compatibility:

- `fastcore.backends.omicverse_*` remain as thin import wrappers.
- `fastcore.vendor.omicverse_gpl` remains as a thin import wrapper.
- `scsp_agent_sop.omicverse_facilities` remains as a thin import wrapper.

This keeps old imports stable while making the OmicVerse license and dependency
boundary explicit. The default Fast environment still does not install the
external OmicVerse package.
