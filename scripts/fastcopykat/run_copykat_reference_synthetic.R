#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: run_copykat_reference_synthetic.R <counts.tsv> <output-dir> <sample-name> [n.cores]", call. = FALSE)
}

counts_path <- args[[1]]
output_dir <- args[[2]]
sample_name <- args[[3]]
n_cores <- if (length(args) >= 4) as.integer(args[[4]]) else 1L

if (dir.exists("/tmp/Rlib-fastcopykat")) {
  .libPaths(c("/tmp/Rlib-fastcopykat", .libPaths()))
}

suppressPackageStartupMessages(library(copykat))

counts <- read.table(counts_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

old <- getwd()
setwd(output_dir)
on.exit(setwd(old), add = TRUE)

set.seed(20260615)
res <- copykat(
  rawmat = counts,
  id.type = "S",
  ngene.chr = 5,
  min.gene.per.cell = 200,
  LOW.DR = 0.05,
  UP.DR = 0.1,
  win.size = 25,
  KS.cut = 0.1,
  sam.name = sample_name,
  distance = "euclidean",
  output.seg = "FALSE",
  plot.genes = "FALSE",
  genome = "hg20",
  n.cores = n_cores
)

cat("result names:", paste(names(res), collapse = ", "), "\n")
cat("prediction dim:", paste(dim(res$prediction), collapse = " x "), "\n")

