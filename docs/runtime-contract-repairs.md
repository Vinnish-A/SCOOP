# Runtime contract repairs

Cross-project failures exposed five shared boundaries. Repairs live at those
boundaries rather than accumulating project-specific fallback rules.

| Boundary | Failure class | Runtime behavior |
| --- | --- | --- |
| Capability versus execution | Dependency discovery succeeded while executable/API or optional algorithm dependency was missing | Doctor identifies discovered adapters/packages and explicitly says execution was not tested; installer imports required RCTD/cNMF APIs; representative jobs validate actual computation |
| Measurement versus library | Whole-transcriptome defaults applied to imaging panels; lineage subsets refit original-library doublets | Explicit assay routing; exact-count QC inheritance; original library IDs control parallel Scrublet jobs |
| File versus object | Backed reads eagerly loaded layers or hit incompatible sparse APIs; export lost usable graph metadata | Compatible AnnData/SciPy bounds and real backed-index smoke; selective read-only CSR/CSC/dense access, metadata-only consumers, preserved descriptors and roundtrip tests |
| Model versus scheduler | Serial independent jobs, hidden train parameters, warmup beyond the budget, nominal NMF seeds with identical starts | Resource-bounded checkpointed jobs; explicit native train/plan options and full histories; independently initialized NMF and actual consensus; genuine full cNMF remains separately named |
| Evidence versus decision | All non-target identities treated as normal reference; candidate-only samples influenced normal scaling | Reviewed eligible reference allowlists, explicit reference errors and reference-only scaling; CNV remains evidence rather than final malignancy classification |

Application content and dependencies have separate identities. Application/Skill
updates reuse locked dependency environments. A genuine dependency change creates
one new dependency identity; immutable old releases remain rollback targets.
The CLI and Python/API launchers bind to the same application and environment.

## Validation scope

Regression tests cover selective reads and order/duplicates, exact-count inherited
QC, real serial/parallel Scrublet equality and input-sensitive resumption, scVI
option precedence/warmup, blocked CNV equivalence, reference-scale isolation,
chromosome-arm normalization and retained graph descriptors after H5AD roundtrip.
Release checks additionally use actual managed CLI/backend jobs, dependency import
smoke, manifest verification and repository tests.

Historical project adapter mistakes and unsupported optional methods are not
retroactively classified as failures of the statistical algorithms. A recorded
bypass is neither native stage completion nor biological approval.

The existing AnnData 0.11 implementation calls the sparse matrix method
`_validate_indices`; SciPy 1.17 moved it to a module-level function. The bound
`scipy<1.17` keeps that existing AnnData generation functional without monkey
patches. The selective reader separately solves eager layer loading; it is not
a substitute for restoring supported dependency combinations. Sources:
[SciPy 1.16 indexing](https://github.com/scipy/scipy/blob/v1.16.3/scipy/sparse/_index.py),
[SciPy 1.17 indexing](https://github.com/scipy/scipy/blob/v1.17.0/scipy/sparse/_index.py).
