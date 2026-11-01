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

免疫、内皮、成纤维、周细胞可作为正常参考候选，但须审阅完整身份、混合来源、技术质量及组织背景；旧标签或 CNV-low 不能直接赋予参考资格。

### 2.4 实际参考、量纲与 checkpoint

同时记录 requested 与 effective reference UID/数量、每 sample 的两套分母、
minimum-reference 参数、被排除参考及原因，以及实际参与 gene selection、
baseline scaling 和 threshold 的 sample scope。外层传入 allowlist 不等于后端
全部采用；少量参考的 sample 可能被默认 minimum 排除，甚至改变选基因范围。
检查后端实际 selection/返回结构，保留未采用参考在完整 query 中的原位置。
说明 baseline 按细胞还是按 sample 加权，不能用传入参考总数替代实际分母。

Checkpoint identity 绑定 count content/UID/gene-order、sample key、reference
label 与实际 UID hash、minimum、实际选基因 sample、gene metadata/build/arm
政策、后端版本/源码和实际参数。配置对象之外的有效参数同样必须绑定。
输入或基线改变时新建 receipt；旧连续结果保留，不能冒称同一 fit。

逐UID桥接成功仍须区分旧/新counts内容及其共同基因比较范围。旧输入版本的
连续profile可作注明版本的支持，不能改称当前原生counts重算，也不能向同群
未覆盖细胞复制。candidate与reference角色可能相交，不能用candidate=False
替代实际reference成员清单。新QC影响多少patient×mode应按真实参考差集核对，
不能据此整体推翻其它有效fit。已有CNV证据的审阅状态与最终恶性assignment
分开记录：支持不足时保留未决，而非一概写成未检查或直接宣告正常/恶性。

区分完整 native depth 与已声明 gene projection 后的 row sum。发现差额先对
省略原列稀疏值和选择政策核验；既不直接宣称 UID 错位，也不放宽检查掩盖
未知差额。分别报告 input columns、genomic common genes、expression-selected
genes、final feature positions、unique symbols 和重复位置政策。不要静默去重
后沿用旧方法名或量纲。

连续 window 直接差异要求 chromosome arm 与 ordered gene list 都一致；相同
window 名称不足以对齐。另行报告 raw-arm 与 thresholded/trimmed-arm 比较、
每 arm 基因覆盖和低覆盖限制；未测 arm 保持 NA。阈值摘要对参考敏感时不能
据此翻转正常/恶性身份。连续 RNA-CNV 敏感性参考资格不自动扩展为空间参考、
最终正常身份或 DNA 证据。分别记录输入准备、算法、导出和核验时间/资源。

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
  triggered_by: adequate reviewed joint evidence remains conflicting or mixed; weak CNV alone is insufficient
not_applicable:
  used_for: immune, stromal, endothelial, mural cells unless tumor lineage is expected
```

Missing evidence is `evidence_pending` workflow progress: actively retrieve or
generate it, document inaccessible evidence and do not submit a final identity or
malignancy decision. This is not a new runtime annotation-schema enum.

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


### Genomic axes, reference execution and algorithm semantics

CNV needs an unambiguous positional axis. Resolve duplicated symbols with verified native stable gene IDs and a build-compatible canonical position table; neither taking the first symbol column nor summing different unknown loci is defensible. Preserve unlocated native feature IDs/counts and the crosswalk. Keep original annotation-version uncertainty explicit when matching versionless stable IDs.

Save requested and actually effective reference UIDs/role labels per sample. A wrapper status contract must accept the status emitted by the real annotation committer; do not rewrite source statuses to manufacture eligibility. The immutable managed release may lag a tested canonical fix, in which case use a recorded validated project adapter. A pooled T/myeloid fit can weight role-specific factors differently from a single merged role, so do not silently change the roles.

Inspect the executed engine: FastCNV can center log2 counts on its expression-selected genomic projection rather than full-library CP10k. Export native depth, genomic counts, unique-selected counts and repeated selected-position counts separately. Ordered selected positions may contain intentional upstream-compatible repeats; they are not the same number of independent genes. Window names alone do not establish comparable ordered genes between legacy and current fits. RNA-CNV evidence and reference sensitivity remain distinct from normal/malignant assignment.
