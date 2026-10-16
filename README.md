# SCOOP

*One scoop is as good as a feast.*

SCOOP means **Single Cell Omics Operating Protocol**. It is a small, auditable workflow for common single-cell and spatial omics projects.

SCOOP follows three principles:

1. **Biology first** — judge identity and malignancy from convergent evidence, not one classifier or CNV threshold.
2. **Measurements first** — preserve patient units, denominators and endpoints; computational success does not establish biology.
3. **Practical execution** — bypass unusable interfaces with a recorded real API or project adapter; resume valid work instead of rebuilding the framework.

Read the [first principles](skill/scoop/references/docs/design-principles.md) for the working contract.

## Workflow

SCOOP is organized as a small set of modules:

```text
qc
core
programs
annotation
spatial
ccc
differential-expression
```

Each module has a clear job.

* `qc`: QC, doublet detection, ambient or suspicious-cell flags.
* `core`: normalization, HVG, PCA, batch correction, graph, UMAP, clustering.
* `programs`: FastCNMF gene program discovery.
* `annotation`: marker, program, CNV, and skill-based evidence for annotation.
* `spatial`: spatial reference and deconvolution when needed.
* `ccc`: candidate cell-cell communication.
* `differential-expression`: pseudobulk condition-level differential expression.

SCOOP avoids adding trajectory, velocity, drug response, foundation-model zoo, or tool-zoo analysis as default modules.

## Fast engines

SCOOP uses a small set of deterministic engines:

```text
FastCore
FastDE
FastCNMF
FastCNVpy
```

They are compute engines.
They do not decide biology by themselves.

Biological labels are adjudicated from tissue-matched atlases, primary literature
and actual expression/sample evidence, then validated and committed deterministically.

## Annotation

SCOOP annotation is evidence-based.

A label is not accepted just because one method says so.
The annotation step collects evidence from:

* marker genes;
* NMF programs;
* curated marker skills;
* reference evidence when available;
* FastCNVpy for tumor samples.

The Python annotation scripts export evidence, prepare a structured decision
template, and commit only validated decisions into H5AD.

Subagent or human edits are allowed, but they must go through the same schema
and validation rules.

## Quick start

```bash
python scripts/install_scoop.py

~/.scoop/bin/scoop init \
  --project <project> \
  --run-id my-run \
  --input-h5ad <project>/data/input.h5ad

~/.scoop/bin/scoop qc \
  --config <project>/runs/my-run/config/run.yaml

~/.scoop/bin/scoop core \
  --config <project>/runs/my-run/config/run.yaml

~/.scoop/bin/scoop programs \
  --config <project>/runs/my-run/config/run.yaml

~/.scoop/bin/scoop annotate evidence \
  --config <project>/runs/my-run/config/run.yaml

~/.scoop/bin/scoop annotate review \
  --config <project>/runs/my-run/config/run.yaml

~/.scoop/bin/scoop annotate commit \
  --config <project>/runs/my-run/config/run.yaml \
  --decisions <project>/runs/my-run/annotation/decisions/annotation_decision_template.json
```

For tumor samples, run FastCNVpy evidence before committing tumor labels:

```bash
~/.scoop/bin/scoop annotate tumor-cnv \
  --config <project>/runs/my-run/config/run.yaml \
  --gene-metadata <project>/data/external/references/gene_metadata.tsv
```

## Repository layout

```text
configs/        default configuration and schemas
docs/           design notes and module documentation
markers/        curated marker skills
scripts/        release helpers and explicit benchmark utilities
skill/          the SCOOP Skill, references, schemas, and report assets
src/            SCOOP and Fast engine code
tests/          regression and validation tests
```

Runtime dependencies have one source of truth: `pyproject.toml`. `uv.lock` is
generated from it and must not be edited manually. The installer creates one
unified environment containing core, atlas, OmicVerse, spatial and CCC
capabilities.
