# Changelog

## Unreleased

- Restore the AnnData/SciPy sparse-index compatibility bound and test real backed slicing during installation.
- Share selective read-only H5AD access across metadata/benchmark readers and add exact-count QC inheritance for retained lineage subsets.
- Parallelize capture-library Scrublet with resource-bounded processes, count/version-bound checkpoints and per-job progress; preserve existing RNA QC and doublet decisions.
- Expose native scVI training/plan options and CLI resource overrides, adapt default KL warmup to the training budget, save training histories and use the shared parallel Leiden sweep for atlas graphs.
- Restrict CNV normal references to eligible reviewed identities, normalize chromosome-arm formats, exclude reference-free samples from reference scaling, and use sparse/blockwise tumor scoring with quadratic clustering opt-in.
- Include official cNMF in the unified environment; report adapter discovery separately from executed validation. Use atomic H5AD writes and external self-file fingerprints; preserve graph/embedding descriptors during delivery and condense the Skill/atlas guidance around measurement, library, execution and review boundaries.

- Separate immutable application releases from dependency environments, so code/Skill updates reuse locked dependencies and rollback restores the corresponding application; a managed Python launcher exposes the same release to API scripts.
- Add a distinct spatial QC/core route with permissive panel-specific defaults, control/coordinate/geometry metrics and preserved counts; retain ordinary scRNA preprocessing and record narrowly scoped LOESS recovery.
- Use actual concurrent resolution/seed Leiden sweeps across managed core backends, with resource-aware workers and a faster two-iteration default; report single-seed stability as unassessed.
- Report the native programme backend as sklearn multistart NMF, use independent random starts and actual consensus components, and retain convergence diagnostics instead of claiming identical deterministic seeds establish stability; expensive OmicVerse validation is explicit.
- Include pinned RCTD-py and call its API directly, preserving full-count exposures and distinguishing sparse multi/doublet fits from their initial full weights. Preserve QC review flags and dry-run/skipped states in stage receipts.

- Allow donor-level abundance and Ro/e comparisons without fitting MIL, so exploratory cell-state screens retain patient-level tests without unnecessary classifiers.

- Added managed donor-level abundance comparisons with exact label permutations, bootstrap intervals, Ro/e, and patient-held-out binary attention MIL with resumable label-permutation validation; composition-only models are explicitly distinguished from expression-based scSurvival.
- Fixed annotation evidence JSON export after H5AD reload and preserved array-backed marker/program names.

- Added opt-in NetworKit ParallelLeiden sweeps with concurrent seed/resolution jobs, OpenMP threads per job, progress timing and graph-bound restart checkpoints; the default Leiden backend is unchanged.
- Added isolated managed-runtime installation without changing the global Skill link.
