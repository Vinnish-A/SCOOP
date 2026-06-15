#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: export_copykat_fixture.R <copykat-source-dir> <output-dir>", call. = FALSE)
}

copykat_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(copykat_dir, "data", "exp.rawdata.rda"))
load(file.path(copykat_dir, "data", "sysdata.rda"))

write.table(exp.rawdata, file.path(output_dir, "copykat_exp_rawdata.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(full.anno, file.path(output_dir, "copykat_full_anno_hg20.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(DNA.hg20, file.path(output_dir, "copykat_DNA_hg20_bins.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
