# Development review and integration record (2026-10-31)

Status: this records a review of SCOOP's development through 2026-10-31, the reconciliation of every branch and
detached tree, and the root-cause repairs that were integrated into `main`. The evidence comes from production use:
the SPI1 glioblastoma atlas (4.2M cells, 182 Visium and 26 CosMx/Xenium sections; issue log
`SPI1/docs/SCOOP_INTERFACE_ISSUES.md`) and the BRCA project (`BRCA/analysis/RUNTIME_DECISIONS.md`). It also uses four
independent code reviews and their repro scripts. Repairs are described as source and test state. They are not a
deployed release and not a scientific revalidation of earlier results.

## 1. Development so far

| Period | Main work |
| --- | --- |
| 2026-06 | Package scaffold: FastCNMF, FastCNVpy, FastDE, the annotation decision layer, FastCore with OmicVerse backends, Harmony 2 core batch correction |
| 2026-08 | SCOOP refactored into an executable Skill; unified SCOOP/OmicVerse runtime with immutable releases and rollback |
| 2026-10-12/16 | Biology-first guidance; shared selective H5AD I/O, parallel Leiden and Scrublet, scVI training options, CNV reference eligibility, atomic writes |
| 2026-10-30 | English Skill edition; lineage-first malignancy guidance (`22fd79b`); evidence/results/recovery design (`docs/08`, `docs/09`) |
| 2026-10-31 | Phase G (`74457eb`, `4ad7f16`, `a5a7601`, `5d67ac5`): measurement mask, mask-aware HVG, typed obs keys, content-bound evidence, declared malignancy, spatial reference producer. Then this review round |

## 2. Branch and tree reconciliation

| Ref or tree | Content | Disposition |
| --- | --- | --- |
| `origin/main` (`bb4dd08`) | Published state before 2026-10-30 evening | Base of everything below |
| local `main` | 4 unpushed commits (`22fd79b`…`a5a7601`) plus a large uncommitted tree: the owner's 2026-10-30 evidence/results/recovery work and the Phase G remainder built on it | Committed as `3e9f290` (owner snapshot, rebuilt on HEAD from the saved owner diff; its tree tested 175/175) and `5d67ac5` (Phase G remainder, 217/217), followed by this round's commits |
| `refs/pull/1/head` (`0044a40`) | Copilot "Initial plan": an empty commit on the initial commit | Nothing to integrate |
| BRCA detached worktree `BRCA/analysis/scoop_runtime_source` at `75259ef` | JSON-safe evidence after an H5AD roundtrip, with its test; an installer flag that skips the global Skill link; RCTD-py 0.3.7; SimpleITK and OpenCV | The first two are on `main` since `1af1d92` (`--no-skill-link`). RCTD is superseded by the locked 0.4.0. SimpleITK and OpenCV serve BRCA's own H&E registration scripts and are not a SCOOP stage, so they stay project-local. The worktree is left in place for BRCA's runtime |
| Installed homes (`~/.scoop`, `.scoop-brca`, `.scoop-ra-parallel`, `.scoop-parallel-validation`) | Built releases, not branches | Unchanged; a new release must be built and activated through `install_scoop.py` |

## 3. What was poorly implemented: root causes

Each theme covers many logged symptoms. The IDs are SPI1 issue-log rows.

1. **Measurement was not modelled.** "Unmeasured" became 0, gene identity was a display string, and normalisers took
   library size from whatever matrix they were given (S034, S044, F02–F08, DC13, MV01–03). Phase G added the
   per-sample mask. Consumers outside QC, core, marker DEG, CNV and spatial still ignore it; they now record
   `not_honoured` instead of staying silent.
2. **Evidence was bound by name, not content.** Commit checked hashes but not the partition, DEG tables were attached
   by prefix, and a missing bundle file made commit rebuild the bundle and then "bind" the partition to itself (S014,
   C1-05, C2-02, FC-06). Archiving re-serialised bundles, so a restored bundle could never match its hash.
3. **Review state had many sources with permissive defaults.** The spatial builder marked every component
   `reviewed`; the adapter hard-coded it; a TSV could override caution; a pilot `gate_override` was ignored (F12/F15,
   S015). The same pattern appeared in annotation, where status words like `not_applicable` slipped past the
   malignancy gate.
4. **There was no stage I/O contract.** Stages found inputs by probing filenames, wrote sidecars that nothing read
   back, bound the whole run.yaml, so any edit made everything stale, and bound no data input when the command
   returned nothing (S022, S023, DC1, C1-08).
5. **Typed keys were partial.** The batch for core Harmony, programme Harmony and scVI each came from a different
   fallback, including the biology-as-batch error, `sample_id` (B03, S010, R10).
6. **Backend contracts were not enforced.** Large inputs reached backends that write no `log1p_norm` (S007); Harmony
   was reported as applied when it did nothing; the legacy core ignored `neighbors_n_pcs`.
7. **Numerical engines produced wrong numbers silently.**
   - FastDE's IRLS used the unfloored working response, so a gene that is zero in one group drifted to log2FC 178
     with P = 0 (S018).
   - FastCNVpy repeated 400 genes in its selection, which added a depth-locked per-cell offset (S042, DC11), and it
     reported windowless arms as neutral (DC5).
   - QC re-admitted cells that an earlier QC had excluded (S026).
   - A float32 FastCNMF prepare crashed consensus (MV05).
8. **Caches and outputs had inconsistent identities.** The RCTD fit identity hashed raw options and the whole QC
   block, but not preset values. Labels did not scope the H5AD or the receipts (F08).
9. **Defaults, docs and copies drifted.** Four hand-synced copies of each packaged YAML had diverged: `configs/`
   advertised `programs.method: fastcnmf` while the code ran sklearn NMF. Docs named Torch Harmony while harmonypy
   ran (S009). The two Skill language editions are maintained by hand.

## 4. Repairs in this round

- **Annotation evidence.**
  - Commit resolves the bundle by the sha256 its decisions cite. It searches beside the decisions, in the decision
    output folder, and in byte-for-byte history.
  - When commit has to build a bundle itself, it reports the partition binding as `not_assessed`.
  - `::group` cannot name a sibling unit, and tables without `label_key` bind nothing.
  - Two malignancy-gate bypasses are closed: `not_applicable` on a tumour-term label, and `ambiguous` at high confidence.
  - Technical exclusions need artefact evidence.
  - A re-commit clears stale malignancy and technical fields.
  - DEG tables are keyed per partition, and the `uns` top-marker side channel is gone.
- **Stage I/O.**
  - Sidecars bind their base by path, hash and obs order, and `read_adata` overlays them onto it.
  - Review, commit and delivery read the object the upstream receipt declared.
  - Receipts bind the config sections each stage consumes and the data file it read.
  - Spatial receipts are kept per label.
  - Delivery keeps exactly what the schema names.
- **Core and QC.**
  - Every backend now writes `log1p_norm`, and the output contract is checked.
  - One batch rule: the declared batch wins, candidates are recorded, and expression-level correction requires a
    declared batch.
  - Harmony reports what actually ran.
  - `neighbors_n_pcs` is honoured.
  - Leiden selection is a declared policy, and the sweep labellings are kept.
  - `size_factor_key` normalises a gene subset by full library size.
  - QC keeps prior decisions as `prior_*`.
- **CNV.**
  - One reference minimum, which makes the q75 opt-in reachable.
  - Per-sample reference flags.
  - Engine configuration from the run config.
  - Mask-only reference gating.
  - Deduplicated gene selection.
  - NaN/`unassessed` for windowless arms.
  - BioMart rules for 5-column metadata.
  - Non-integer X is refused as counts.
- **Engines.**
  - FastDE uses DESeq2's working response; log2FC matches R DESeq2 1.46 to within 4e-4, and diverged fits are flagged.
  - FastCNMF keeps dtypes, applies a convergence rule and records its ridge layout.
  - `validate` refuses duplicate IDs and non-integer counts and returns outputs.
  - `ccc` tests only accepted labels.
  - Atlas saves the model before the graph steps, takes `model_kwargs` and `hvg_key`, records the chosen resolution,
    refuses absent covariates, builds log1p markers and registers label-keyed DEG.
- **Spatial.**
  - One review-status derivation, defaulting to `evidence_pending`, with per-component interpretation.
  - One native-exposure rule shared by builder, adapter and harness.
  - Donor-proportional holdouts.
  - A fit identity built from resolved compute settings.
  - Label-bound artifacts.
  - Dry-run checks every input without fitting.
  - Identifiability tables bound to their reference.
- **Configuration.** One packaged source for each YAML, with byte-identical copies enforced by a test. The scVI
  covariates default to none.

Validation:

- The full repository suite passes on the integrated tree (302 tests), run as
  `PYTHONPATH=src <runtime env>/bin/python -m pytest tests`. Each themed commit's staged tree passed the suite on its own.
- Each repair carries a regression test; the engine tests were also confirmed to fail on the earlier code.
- A wheel built from the committed tree contains the new modules and the packaged YAMLs, and its CLI starts.
- A CLI smoke run on synthetic counts ran validate, qc, core, programs, annotate evidence/review/commit and
  deliver. It used `obs-parquet` and kept decisions outside the review folder.
  - Commit found the review's bundle by hash. Deliver read the sidecar the commit receipt declared.
  - `status --verify` reported nothing stale after an unrelated config edit, and only qc, core and validate after
    a QC edit.
  - The run exposed the crash on unscorable programme gene sets, which is now fixed.

Behaviour changes users will see:

- RCTD fits refit once, because the fit identity changed.
- Production spatial runs refuse references without committed review.
- Atlas runs must declare their covariates.
- The CNV threshold needs 20 references.
- Legacy core uses 30 neighbour PCs.

## 5. Still open, ranked by consequence

1. **Whole-object I/O.** Stages still read and write full H5ADs (S023, S028, S032, G01, G08, MV04). Sidecars now
   bind correctly, but there is no streaming row/column data layer.
2. **Measurement in model fitting.** There is no masked likelihood for scVI (F08, DV2-03), no mask in pseudobulk DE
   or FastCNMF (DV2-10), and no command that produces the measurement mask.
3. **CNV engine.** No fitted-reference object for fit and apply (DV3, DV6, DV8). No focal resolution (DC6). No per-arm
   top-up, which would need windows built from the final gene set. Threshold scope is not recorded (DV5). Joins use
   symbols rather than stable IDs (S041, DV2).
4. **FastDE.** No covariates (`~ batch + condition`), no stage receipt, no Cook's distance and no independent
   filtering; dispersion outlier spread uses a variance where DESeq2 uses MAD².
5. **Annotation outcomes.** Still to do:
   - pending units (C1-02);
   - broad-scope fine status (C1-03, owner decision);
   - citable `register-result` artefacts (C1-04);
   - modality-limited evidence (E03);
   - structured unit lineage (C2-01);
   - a reviewed-delivery route for the atlas pipeline.
6. **Spatial.** Imaging doublet chunking (F14), a merged-package command (F13), a content-addressed reference cache
   (M02) and a per-section loop (S005). Bulk deconvolution, survival models, IMC, ATAC and CITE are project adapters,
   not stages.
7. **Multi-key declared batch.** `run.obs_keys.batch` is a single column, and FastCNMF HVGs ignore batch (S012, G06).
8. **Skill editions.** The original-language edition lags the English one; for example, `stages/spatial.md` still
   names RCTD 0.3.9. The editions should be generated from one source or checked against each other.
