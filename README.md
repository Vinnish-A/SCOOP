# SCOOP

*One scoop is as good as a feast.*

SCOOP means **Single Cell Omics Operating Protocol**. It is a small, auditable workflow for common single-cell and spatial omics projects.

SCOOP follows three rules:

1. **Fast** — use lightweight deterministic engines where possible.
2. **Robust** — do not trust a single marker, model, or reference.
3. **Simple** — keep one default method; keep fallbacks explicit.

In other words: do not boil the ocean. Take one validated scoop.

## Workflow

SCOOP is organized as a small set of modules:

```text
01_qc
02_core
03_programs
04_annotation
05_spatial
06_ccc
07_de
```

Each module has a clear job.

* `01_qc`: QC, doublet detection, ambient or suspicious-cell flags.
* `02_core`: normalization, HVG, PCA, batch correction, graph, UMAP, clustering.
* `03_programs`: FastCNMF gene program discovery.
* `04_annotation`: marker, program, CNV, and skill-based evidence for annotation.
* `05_spatial`: spatial reference and deconvolution when needed.
* `06_ccc`: candidate cell-cell communication.
* `07_de`: pseudobulk condition-level differential expression.

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

Biological labels are proposed by a subagent or analyst using world knowledge,
then checked by the annotation schema and validator before they are committed.

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
pip install -e .

mkdir -p .scoop_local/runs .scoop_local/data
cp -r templates/run .scoop_local/runs/my_run

python scripts/01_qc_scrublet.py \
  --config .scoop_local/runs/my_run/config/run.yaml

python scripts/02_core_analysis.py \
  --config .scoop_local/runs/my_run/config/run.yaml

python scripts/03_fast_consensus_nmf.py \
  --config .scoop_local/runs/my_run/config/run.yaml

python scripts/04a_annotation_evidence.py \
  --config .scoop_local/runs/my_run/config/run.yaml

python scripts/04c_annotation_export_for_agent.py \
  --config .scoop_local/runs/my_run/config/run.yaml

python scripts/04d_annotation_commit.py \
  --config .scoop_local/runs/my_run/config/run.yaml \
  --decisions .scoop_local/runs/my_run/04_annotation/decisions/annotation_decision_template.json
```

For tumor samples, run FastCNVpy evidence before committing tumor labels:

```bash
python scripts/04b_tumor_fastcnvpy.py \
  --config .scoop_local/runs/my_run/config/run.yaml \
  --gene-metadata .scoop_local/data/external/references/gene_metadata.tsv
```

## Repository layout

```text
configs/        default configuration and schemas
docs/           design notes and module documentation
markers/        curated marker skills
scripts/        executable workflow entry points
src/            SCOOP and Fast engine code
tests/          regression and validation tests
```
