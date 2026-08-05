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
broad cell types -> candidate malignant compartment -> FastCNVpy -> malignant_status
```

候选 malignant compartments 通常包括：

```text
epithelial/parenchymal/glial/neural crest-derived tumor cells
```

免疫、内皮、成纤维、周细胞通常作为 CNV-neutral reference，但需要排除 doublet 和低质量细胞。

### 2.2 当正常与肿瘤分得开

```text
marker tumor-like + CNV abnormal -> malignant-high
marker normal-like + CNV neutral -> normal-high
```

### 2.3 当正常与肿瘤分不开

```text
same broad epithelial/parenchymal cluster
  -> CNV abnormal arm-level pattern: malignant/probable malignant
  -> CNV neutral and normal-reference matched: normal parenchymal
  -> CNV weak/conflicting: ambiguous, keep low confidence
```

---

## 3. malignant evidence hierarchy

| evidence level | evidence type | interpretation |
|---|---|---|
| A | patient-specific CNV / CNA from FastCNVpy, inferCNV, CopyKAT, Numbat or equivalent | strongest scRNA-level evidence |
| A | matched WES/WGS/panel/scDNA/FISH/IHC mutation or amplification evidence | strongest external evidence |
| B | patient-specific tumor cluster plus poor normal-reference mapping | supportive |
| B | cancer-specific malignant program plus tumor-region enrichment | supportive |
| C | only canonical epithelial/parenchymal markers | not enough |
| C | only proliferation or stress | not enough |

---

## 4. malignancy_status labels

```yaml
malignant-high:
  required: CNV abnormal or tumor SNV positive
malignant-probable:
  required: weak/partial CNV or tumor-specific cluster plus poor normal-reference mapping
normal:
  required: CNV neutral and normal-reference matched
ambiguous:
  triggered_by: low quality, doublet risk, weak CNV, normal and tumor programs mixed
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
Do not call functional states before malignant_status is assigned.
Do not interpret CNV state as functional state; CNV separates malignant/normal and clones, not cNMF programs.
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
