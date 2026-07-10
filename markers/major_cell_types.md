# 大群细胞类型注释 marker 词典

版本：v1.0-refined  
定位：用于全体细胞的 broad annotation。该文件只做大类，不做精细亚型。

---

## 1. 大群注释顺序

```text
1. immune vs non-immune
2. epithelial / parenchymal vs stromal vs endothelial vs mural
3. immune 内部：T/NK/B/plasma/myeloid/mast/neutrophil
4. tumor sample: 对可疑实质/上皮/胶质群运行 FastCNVpy
5. 输出 broad_type 和需要进入哪个精细 skill
```

---

## 2. 大群 marker 表

| broad_type | positive markers | exclusion markers | next skill |
|---|---|---|---|
| `T` | `CD3D, CD3E, CD3G, TRAC, TRBC1, TRBC2, CD2, CD7` | `MS4A1, CD79A, LST1, C1QA, EPCAM, COL1A1` | `T.md` |
| `NK` | `NKG7, GNLY, PRF1, GZMB, KLRD1, KLRF1, FCGR3A, NCR1` | `CD3D, CD3E, TRAC` for pure NK | project NK extension |
| `B` | `MS4A1, CD79A, CD79B, CD74, HLA-DRA, CD19, CD22` | `CD3D, LST1, EPCAM` | B extension |
| `Plasma` | `MZB1, XBP1, JCHAIN, IGKC, IGHG1, IGHA1, SDC1` | strong T/myeloid markers | B extension |
| `Myeloid` | `LYZ, LST1, TYROBP, FCER1G, AIF1, CST3, CTSS` | `CD3D, MS4A1, EPCAM, COL1A1` | `myeloid.md` |
| `Monocyte` | `FCN1, VCAN, S100A8, S100A9, CCR2, LILRB1` | strong epithelial/stromal | `myeloid.md` |
| `Macrophage` | `C1QA, C1QB, C1QC, APOE, MERTK, MSR1, CD68` | `CLEC9A/XCR1` cDC1, `CD1C/FCER1A` cDC2 if dominant | `myeloid.md` |
| `DC` | `FCER1A, CD1C, CLEC10A, CLEC9A, XCR1, LAMP3, CCR7` | macrophage genes if stronger | `myeloid.md` |
| `Mast` | `KIT, TPSAB1, TPSB2, CPA3, MS4A2` | basophil/neutrophil markers | `myeloid.md` |
| `Neutrophil` | `S100A8, S100A9, CSF3R, CXCR2, FCGR3B, MPO, ELANE` | macrophage/DC markers if dominant | `myeloid.md` |
| `Epithelial` | `EPCAM, TACSTD2, KRT8, KRT18, KRT19, CDH1, CLDN3, CLDN4` | `PTPRC, COL1A1, PECAM1, RGS5` | `tumor-general.md` or tissue skill |
| `MalignantEpi` | epithelial markers plus CNV/SNV/tumor reference evidence | CNV-neutral normal reference | cancer-specific tumor skill |
| `Fibroblast` | `COL1A1, COL1A2, COL3A1, DCN, LUM, COL6A1, DPT, PDGFRA` | `PTPRC, EPCAM, PECAM1, RGS5, MYH11` | `fibroblast.md` |
| `CAF` | fibroblast markers plus `FAP, PDPN, POSTN, MMP11, COL11A1, CTHRC1, LRRC15` | immune/endothelial/mural | `fibroblast.md` |
| `Endothelial` | `PECAM1, VWF, KDR, FLT1, CLDN5, EMCN, PLVAP, ESAM` | `RGS5, ACTA2, COL1A1, PTPRC` | endothelial extension |
| `LymphaticEndo` | `PROX1, PDPN, LYVE1, FLT4, CCL21` | macrophage `LYVE1` alone is insufficient | endothelial extension |
| `Pericyte` | `RGS5, CSPG4, PDGFRB, MCAM, NOTCH3` | endothelial `PECAM1/VWF`, fibroblast `DCN/LUM` | mural extension |
| `SMC` | `ACTA2, TAGLN, MYH11, CNN1, MYL9` | pericyte/fibroblast context | mural extension |
| `Erythroid` | `HBB, HBA1, HBA2, ALAS2, SLC4A1` | most tissue markers | remove/review unless expected |
| `Platelet` | `PPBP, PF4, GP9, ITGA2B, NRGN` | nucleated cell markers | remove/review unless expected |

---

## 3. 冲突规则

### 3.1 doublet 优先排除

```text
PTPRC + EPCAM/KRT high -> immune-epithelial doublet risk
PTPRC + COL1A1/DCN high -> immune-stromal doublet risk
PECAM1/VWF + RGS5/MYH11 high -> endothelial-mural doublet or vascular-associated capture
CD3D/TRAC + NKG7/GNLY high -> T/NK-like T or NK-T ambiguity; check TRAC/CD3 strength
HBB/HBA/PPBP/PF4 high across many lineages -> ambient blood/platelet RNA
```

### 3.2 marker 不可单独使用

```text
EPCAM alone does not prove malignant epithelial.
COL1A1 alone does not prove fibroblast; it may be ambient or EMT tumor.
LYZ alone does not prove myeloid subtype.
CD74/HLA-DRA alone does not prove B cell, DC, macrophage or apCAF.
MKI67/TOP2A only indicate cycling.
```

---

## 4. 输出

```yaml
broad_annotation:
  broad_type: T | NK | B | Plasma | Myeloid | Epithelial | Fibroblast | Endothelial | Mural | Erythroid | Platelet | Unknown
  positive_evidence: []
  exclusion_evidence: []
  doublet_risk: none | low | medium | high
  next_skill: <file>
  confidence: high | medium | low
```
