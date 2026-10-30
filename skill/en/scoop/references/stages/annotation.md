# Annotation evidence


Inputs:

- Credible clustering and conventional DEG with tissue marker evidence for broad review.
- For formal fine-subtype review, the accepted biological parent and all evidence
  required by the [annotation workflow](../markers/annotation-workflow.md),
  including batch-handled parent NMF and within-parent clustering/DEG.

Script:

```bash
scoop annotate evidence --config runs/<run_id>/config/run.yaml
```

OmicVerse fast marker option:

```bash
scoop annotate evidence --config runs/<run_id>/config/run.yaml --use-omicverse-cosg
```

The OmicVerse wrapper only exports marker evidence. Delegates may assemble evidence
and propose labels; the main agent personally reads the actual evidence and makes
each final biological decision, which must pass validation before commit. Broad
review does not require NMF; an incomplete fine review remains provisional.

H5AD writes:

- No final labels without the main agent's evidence-based adjudication and validation.

Sidecar outputs:

- `tables/cluster_markers.parquet`.
- `tables/annotation_evidence_template.tsv`.

Review triggers:

- Reference and marker evidence disagree.
- New cell type claim.
- Low-confidence annotation used in downstream conclusions.
