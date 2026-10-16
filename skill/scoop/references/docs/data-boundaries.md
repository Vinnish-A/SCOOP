# Matrix, library and persistence boundaries

## Selective H5AD reads

AnnData backed mode may eagerly load layers, and some AnnData/SciPy combinations
break backed sparse slicing. The managed dependency contract pins SciPy below
1.17 while AnnData remains below 0.12 and tests actual backed slicing at install.
Use the managed reader when selecting a small subset
or inspecting metadata, rather than loading all layers or changing the shared stack:

```python
from scsp_agent_sop.h5ad_io import H5ADReader

with H5ADReader("input.h5ad") as source:
    obs, var = source.obs, source.var
    counts = source.read_rows([10, 2, 10], layer="counts")
```

CSR, CSC and dense matrices are supported. Requests retain order and duplicates;
empty requests are valid. CSR row ranges are read together; CSC reads selected
columns then selects rows. Choose block sizes for the matrix layout and memory.
This is selective access, not a claim that every downstream algorithm is out-of-core.
The reader opens source files read-only and does not alter the dependency runtime.

## Library-level computation and subsets

Doublet simulation belongs to original capture libraries. Configure `keys.sample`
to the capture-library ID, not the donor or a biological compartment. Ordinary RNA
keeps the existing QC/Scrublet rules. Independent libraries run in separate bounded
processes, with one thread each by default, per-job timing and count/ID/parameter/
version-bound checkpoints. Failed work remains inspectable; completed libraries
resume without refitting.

For a lineage extracted from a QC'd library, inherit decisions:

```bash
scoop qc --config run.yaml --input lineage.h5ad --inherit-from parent_qc.h5ad
```

The runtime verifies unique IDs, unchanged gene axis, library membership, retained
parent cells and exact raw counts before copying decisions. It records parent SHA
and inherited scope, without claiming new QC/doublet estimation. A feature-filtered
object is not an exact full-count subset; recover the source gene axis first.

## Persistence

Keep counts and feature-measurement meaning intact. Repair field types at the
source metadata boundary; do not stringify all numeric data or silently drop
unserializable evidence. Use JSON for nested list-of-record parameter payloads and
external tables for large results. Managed QC/core/programme/spatial/atlas/delivery writes replace the output only
after successful serialization. The embedded registry does not claim a file
contains its own hash; the external manifest fingerprints the final object.
Verify a saved/read-back object rather than
repeating model fitting after an export failure.

Pruning retains graph descriptors that point to retained adjacency/distance
matrices and embedding descriptors for retained PCA/UMAP coordinates. Use a
project delivery schema for additional clinical fields and custom graph names.
