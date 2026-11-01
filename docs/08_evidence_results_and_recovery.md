# Evidence, usable results and recovery

Status: implementation approved by the 2026-10-30 SPI1 user request; this document
describes the intended contract. Execution and validation receipts determine what
is actually implemented and deployed.

## Problem

The same word, “complete”, currently covers process exit, saved output, numerical
validation, reviewed identity and supported interpretation. Discovery in doctor,
an H5AD filename, an issue note or a canonical source fix can then be mistaken for
a usable analysis or deployed repair. Project adapters repeat count/QC joins,
reference construction and recovery logic because the normal interface cannot
consume their verified artifacts efficiently.

The current workflow wrapper can mark a no-output command committed, inventories
old files as if they were new stage output, omits input bindings and overwrites
failed receipts. Delivery can silently fall back to evidence/program/QC input.
The RCTD adapter handles only final_use, can drop rare reference types, discards
excluded observations from delivery and cannot resume export after a valid fit.

## Design

Use three existing boundaries rather than a new execution framework:

- A measurement preserves its native axes, count scale, full exposure and source
  selection. Explicit source failures apply before fitting; absence of one QC
  field does not erase separate QC evidence. Missing measurements stay missing.
- A result records computation, validation, interpretation and coverage
  separately. An external project adapter can register the same artifact contract
  with its actual route and checks, without pretending a managed stage ran.
- Delivery exports the declared result scope. Exporting current state remains
  possible explicitly; it does not silently become reviewed annotation. Fit
  checkpoints and export are separate so a serialization repair reuses valid work.

Keep existing stage commands and configs. Add a readable status command and a
bounded result-registration route. Verify declared hashes on request; show stale,
missing, unassessed and limited results separately. Do not require every unrelated
cell to resolve before independent computations can run. Unsupported fine
interpretation does not make a successful exploratory fit disappear.

For large references, accept a hash-bound native-rate profile cache with explicit
review and exposure provenance. Preserve every declared reference component;
report rare/uncertain types instead of silently dropping them. Keep fast and
standard RCTD presets with actual consumed parameters. Upgrade the known older
calibration implementation only through the locked runtime release workflow.

## Implementation and verification

1. Repair CLI discovery, receipt truthfulness/history, result registration and
   status inspection. A successful no-output/dry-run command cannot imply a
   validated result; stale inputs/outputs must be visible.
2. Require explicit partial delivery, retain reviewed unresolved cells, and bind
   the delivered scope to the selected input and outputs.
3. Repair spatial source masks, count/exposure handling, rare-type coverage,
   cached references and fit/export recovery. Keep all observation rows with NA
   for excluded/unfitted weights. Validate against a real SPI1 section.
4. Build a locked candidate, run meaningful regressions and the repository suite,
   test a cold installation, then activate the same content ID locally with a
   rollback target. Do not report source tests as active runtime validation.

Scientific repairs are recorded in SPI1 separately: finish the four omitted
current CNV modes; quantify unresolved annotation coverage; compare deconvolution
under matched reference/feature conditions. Registration and new code do not
retroactively validate those analyses.

No new remote service, framework dependency, automatic biological classifier,
blanket approval gate, Git commit or publication is part of this change.

## Current interfaces

`scoop --help` exposes the semantic commands already routed by the CLI.
`scoop register-result --run-root RUN --manifest result.json` binds actual
project-adapter inputs/outputs and declared checks; `scoop status --run-root RUN
--verify` inspects them without equating file integrity with scientific truth.
The shared contract is documented in the Skill's workflow reference. A registrar
cannot detect omitted dependencies or replace independent scientific review.

Annotation exports all tested input genes, keeps auxiliary specificity scores
separate, and connects registered DEG to decision templates. Missing DEG citation
requires review; cited registered files must remain present/hash-matched at the
commit boundary. `identity_scope` expresses broad, fine or unresolved decisions.
Explicit unresolved labels can be recorded with review_required status, while a
rejected replacement cannot retain accepted eligibility from an earlier result.
Undeclared identity scope defaults to broad. Fine acceptance requires an explicit
fine declaration; legacy accepted/RNA-phenotype flags do not establish subtype
resolution. A scope-only metadata correction can reuse models after exact checks
show their consumed labels, eligibility, measurements and axes did not change.

The spatial API accepts native H5AD counts or a reviewed native-rate profile
cache with its manifest. Source registry joins validate exact row/ID axes. Count
selection, full native exposure, reference component retention and original-row
delivery are separate invariants. A fit-identity checkpoint precedes export and
can recover from a serialization failure. An explicit all_shared-gene sensitivity
supports controlled comparisons; default fitting still selects DE genes.

The release evidence belongs in the consuming project's deployment receipt:
cold locked installation/import smoke, actual backend pilot, regression suite,
manifest hashes, final candidate/production content IDs and rollback target.
Source changes and those checks must not be described as active deployment until
the actual managed links and package/module paths have been verified.
