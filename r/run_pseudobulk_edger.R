#!/usr/bin/env Rscript
# Run edgeR quasi-likelihood pseudobulk DE for one cell-type directory.
#
# Usage:
#   Rscript r/run_pseudobulk_edger.R runs/<run_id>/07_de/pseudobulk/<cell_type> condition ctrl test
#
# Required files in the cell-type directory:
#   counts.tsv: rows=pseudobulk samples, columns=genes, first column pseudobulk_id
#   metadata.tsv: rows=pseudobulk samples, must include pseudobulk_id and condition column

suppressPackageStartupMessages({
  library(edgeR)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript r/run_pseudobulk_edger.R <pseudobulk_dir> <condition_col> <ctrl_group> <test_group>")
}
indir <- args[[1]]
condition_col <- args[[2]]
ctrl <- args[[3]]
test <- args[[4]]
counts_path <- file.path(indir, "counts.tsv")
meta_path <- file.path(indir, "metadata.tsv")
counts <- read.delim(counts_path, check.names=FALSE)
meta <- read.delim(meta_path, check.names=FALSE)
rownames(counts) <- counts$pseudobulk_id
counts$pseudobulk_id <- NULL
rownames(meta) <- meta$pseudobulk_id
common <- intersect(rownames(counts), rownames(meta))
counts <- counts[common, , drop=FALSE]
meta <- meta[common, , drop=FALSE]
meta <- meta[meta[[condition_col]] %in% c(ctrl, test), , drop=FALSE]
counts <- counts[rownames(meta), , drop=FALSE]
meta[[condition_col]] <- factor(meta[[condition_col]], levels=c(ctrl, test))
if (min(table(meta[[condition_col]])) < 2) {
  warning("Fewer than two pseudobulk samples in at least one group; results are exploratory.")
}
y <- DGEList(counts=t(counts), samples=meta)
keep <- filterByExpr(y, group=meta[[condition_col]])
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(as.formula(paste0("~", condition_col)), data=meta)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
qlf <- glmQLFTest(fit, coef=ncol(design))
res <- topTags(qlf, n=Inf)$table
res$gene <- rownames(res)
outdir <- file.path(dirname(dirname(indir)), "contrasts", paste0(test, "_vs_", ctrl), basename(indir))
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
write.table(res, file=file.path(outdir, "de_edgeR.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(design, file=file.path(outdir, "design_matrix.tsv"), sep="\t", quote=FALSE, row.names=TRUE)
