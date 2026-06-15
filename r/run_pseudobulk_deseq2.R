#!/usr/bin/env Rscript
# Run DESeq2 Wald pseudobulk DE for one cell-type directory.

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript r/run_pseudobulk_deseq2.R <pseudobulk_dir> <condition_col> <ctrl_group> <test_group>", call. = FALSE)
}
if (!requireNamespace("DESeq2", quietly=TRUE)) {
  stop("DESeq2 is not installed", call. = FALSE)
}

indir <- args[[1]]
condition_col <- args[[2]]
ctrl <- args[[3]]
test <- args[[4]]
counts <- read.delim(file.path(indir, "counts.tsv"), check.names=FALSE)
meta <- read.delim(file.path(indir, "metadata.tsv"), check.names=FALSE)
rownames(counts) <- counts$pseudobulk_id
counts$pseudobulk_id <- NULL
rownames(meta) <- meta$pseudobulk_id
common <- intersect(rownames(counts), rownames(meta))
counts <- counts[common, , drop=FALSE]
meta <- meta[common, , drop=FALSE]
meta <- meta[meta[[condition_col]] %in% c(ctrl, test), , drop=FALSE]
counts <- counts[rownames(meta), , drop=FALSE]
meta[[condition_col]] <- factor(meta[[condition_col]], levels=c(ctrl, test))

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = t(round(as.matrix(counts))),
  colData = meta,
  design = as.formula(paste0("~", condition_col))
)
dds <- dds[rowSums(DESeq2::counts(dds)) >= 10, ]
dds <- DESeq2::DESeq(dds, quiet=TRUE)
res <- as.data.frame(DESeq2::results(dds, contrast=c(condition_col, test, ctrl)))
res$gene <- rownames(res)
outdir <- file.path(dirname(dirname(indir)), "contrasts", paste0(test, "_vs_", ctrl), basename(indir))
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
write.table(res, file=file.path(outdir, "de_DESeq2.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
