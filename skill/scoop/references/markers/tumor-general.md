# 肿瘤组织通用 malignant / normal 判定 skill

版本：v1.0-refined  
范围：单细胞转录组学肿瘤样本。  
定位：在进入癌种专属 tumor skill 前，先判断细胞是否为 malignant。

---

## 1. 总原则

```text
先证明 malignant，再做癌种内部亚型。
```

在肿瘤组织中，正常实质细胞和肿瘤细胞常共享 marker。例如：

- 乳腺正常 luminal 与 luminal breast cancer 都可表达 `EPCAM/KRT8/KRT18/KRT19/GATA3`。
- CRC 正常 crypt/absorptive/goblet 与肿瘤上皮都可表达肠上皮 marker。
- GBM 正常 astrocyte/OPC/neuron 与 AC-like/OPC-like/NPC-like tumor marker 高度重叠。

因此，marker 只能确定“候选来源”，不能单独确定 malignant。

---

## 2. FastCNVpy 使用位置

### 2.1 大群阶段

```text
broad cell types -> candidate compartment -> CNV/SNV/reference/context review -> malignant_status
```

候选 malignant compartments 通常包括：

```text
epithelial/parenchymal/glial/neural crest-derived tumor cells
```

免疫、内皮、成纤维、周细胞通常作为 CNV-neutral reference，但需要排除 doublet 和低质量细胞。

### 2.2 当正常与肿瘤分得开

```text
tumor-like programme + concordant CNV/context evidence -> reviewed malignant support
coherent normal lineage + normal-atlas and tissue-context support -> normal-supported
CNV neutral alone -> no automatic normal label
```

### 2.3 当正常与肿瘤分不开

```text
same broad epithelial/parenchymal cluster
  -> CNV abnormal arm-level pattern: supporting evidence, adjudicate with identity/reference/context
  -> coherent normal lineage plus normal-reference and tissue-context support: normal parenchymal (CNV neutrality alone is insufficient)
  -> CNV low/diploid + concordant tumor programme, patient-specific structure and tumor context: malignant-probable remains eligible
  -> genuinely conflicting or insufficient multi-source evidence: ambiguous
```

---

## 3. malignant evidence hierarchy

| evidence level | evidence type | interpretation |
|---|---|---|
| B | expression-inferred CNV / CNA from FastCNVpy, inferCNV, CopyKAT, Numbat or equivalent | supporting RNA evidence; sensitive to reference, depth and state, not independent DNA confirmation |
| A | matched WES/WGS/panel/scDNA/FISH/IHC mutation or amplification evidence | strongest external evidence |
| B | patient-specific tumor cluster plus poor normal-reference mapping | supportive |
| B | cancer-specific malignant program plus tumor-region enrichment | supportive |
| C | only canonical epithelial/parenchymal markers | not enough |
| C | only proliferation or stress | not enough |

---

## 4. malignancy_status labels

```yaml
malignant-high:
  required: convergent tumor evidence and high-confidence review; RNA-CNV alone is insufficient and high confidence is not DNA confirmation
malignant-probable:
  required: convergent tumor-context, patient-specific structure and tumor-program/normal-reference evidence; weak or absent CNV does not exclude this label
normal:
  required: coherent normal lineage plus matched-normal-reference and tissue/context support; CNV neutrality alone is insufficient
ambiguous:
  triggered_by: insufficient/conflicting joint evidence, including unresolved QC or mixed-origin risks; weak CNV alone is insufficient
not_applicable:
  used_for: immune, stromal, endothelial, mural cells unless tumor lineage is expected
```

---

## 5. epithelial / parenchymal broad gate

| tissue context | candidate cell markers | cancer skill |
|---|---|---|
| breast | `EPCAM, TACSTD2, KRT8, KRT18, KRT19, CDH1, MUC1` | `tumor-BRCA.md` |
| colorectal | `EPCAM, TACSTD2, KRT8, KRT18, KRT19, KRT20, CDH1, CLDN3, CLDN4` | `tumor-CRC.md` |
| GBM / brain | `GFAP/AQP4/ALDOC` astrocytic, `PDGFRA/OLIG1/OLIG2` OPC-like, `SOX2/SOX4/ASCL1/DCX` NPC-like, `CD44/CHI3L1/VIM` MES-like | `tumor-GBM.md` |

---

## 6. do_not_call

```text
Do not call malignant from EPCAM/KRT/GFAP/PDGFRA alone.
Do not call normal from CNV-negative alone if tumor type can be CNV-low or MSI-like; use SNV/reference/tumor program.
Keep differentiation and functional states available while malignant_status is assessed separately; unresolved malignancy must not erase a supported epithelial identity.
Do not turn RNA CNV thresholds into final malignant/normal labels. Retain continuous CNV/arm evidence separately; inspect reference, low-RNA, cycling and stress effects, and resolve malignancy using the complete biological evidence.
Do not let cycling/ribosomal/stress clusters become tumor subtypes.
```

---

## 7. recommended output

```yaml
tumor_status_record:
  cluster_id: <cluster>
  broad_type: epithelial_or_parenchymal
  cnv_status: abnormal | neutral | weak | not_run
  snv_status: positive | negative | unknown
  normal_reference_mapping: matched | poor | unknown
  malignant_status: malignant-high | malignant-probable | normal | ambiguous
  next_skill: tumor-BRCA.md | tumor-CRC.md | tumor-GBM.md | normal tissue annotation
  confidence: high | medium | low
```
