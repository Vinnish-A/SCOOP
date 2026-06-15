#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: run_r_cnvcalling_reference.R <counts.tsv> <geneMetadata.tsv> <output-dir> [obs.tsv] [referenceVar] [referenceLabel]", call. = FALSE)
}

counts_path <- args[[1]]
gene_metadata_path <- args[[2]]
output_dir <- args[[3]]
obs_path <- if (length(args) >= 4 && nzchar(args[[4]])) args[[4]] else NA
reference_var <- if (length(args) >= 5 && nzchar(args[[5]])) args[[5]] else NA
reference_label <- if (length(args) >= 6 && nzchar(args[[6]])) strsplit(args[[6]], ",", fixed = TRUE)[[1]] else NA

threshold_percentile <- 0.01
window_size <- 150
window_step <- 10
top_n_genes <- 7000

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rawCounts <- as.matrix(read.table(counts_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
geneMetadata <- read.table(gene_metadata_path, sep = "\t", header = TRUE, check.names = FALSE, quote = "", comment.char = "")
obs <- NULL
if (!is.na(obs_path)) {
  obs <- read.table(obs_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, quote = "", comment.char = "")
}

scale_on_reference_label <- TRUE
if (is.na(reference_var) || is.na(reference_label[1]) || is.null(obs)) {
  scale_on_reference_label <- FALSE
  referenceCells <- colnames(rawCounts)
} else {
  if (length(reference_label) == 1) {
    referenceCells <- rownames(obs)[which(obs[[reference_var]] == reference_label)]
    if (length(referenceCells) == 0) {
      scale_on_reference_label <- FALSE
      referenceCells <- colnames(rawCounts)
    }
  } else {
    referenceCells <- list()
    for (label in reference_label) {
      cells <- rownames(obs)[which(obs[[reference_var]] == label)]
      if (length(cells) >= 5) referenceCells[[label]] <- cells
    }
    if (length(referenceCells) == 0) {
      scale_on_reference_label <- FALSE
      referenceCells <- colnames(rawCounts)
    }
  }
}

geneMetadata <- geneMetadata[which(geneMetadata$gene_biotype %in% c("protein_coding", "lncRNA") &
                                     geneMetadata$chromosome_name %in% c(1:22, "X") &
                                     geneMetadata$hgnc_symbol != ""), ]
geneMetadata$chromosome_num <- geneMetadata$chromosome_name
geneMetadata$chromosome_num[which(geneMetadata$chromosome_num == "X")] <- 23
geneMetadata$chromosome_num <- as.numeric(geneMetadata$chromosome_num)
geneMetadata2 <- unique(geneMetadata[, c("hgnc_symbol", "chromosome_num", "start_position", "end_position", "chr_arm")])

funTrim <- function(normcounts, lo = -3, up = 3) {
  t(apply(normcounts, 1, function(z) {
    z[which(z < lo)] <- lo
    z[which(z > up)] <- up
    z
  }))
}

commonGenes <- intersect(rownames(rawCounts), geneMetadata2$hgnc_symbol)
rawCounts <- rawCounts[commonGenes, , drop = FALSE]

if (scale_on_reference_label) {
  if (is.list(referenceCells)) {
    allRef <- unlist(referenceCells)
    averageExpression <- rowMeans(rawCounts[, allRef, drop = FALSE])
  } else {
    averageExpression <- rowMeans(rawCounts[, referenceCells, drop = FALSE])
  }
} else {
  averageExpression <- rowMeans(rawCounts)
}

if (length(commonGenes) < top_n_genes) top_n_genes <- length(commonGenes)
topExprGenes <- commonGenes[order(averageExpression, decreasing = TRUE)[1:top_n_genes]]

topExprGenes_metadata <- geneMetadata2[geneMetadata2$hgnc_symbol %in% topExprGenes, ]
topExprGenes_metadata$chr_arm_full <- paste0(topExprGenes_metadata$chromosome_num, topExprGenes_metadata$chr_arm)

genes_by_arm <- split(topExprGenes_metadata$hgnc_symbol, topExprGenes_metadata$chr_arm_full)

for (arm in unique(geneMetadata2$chr_arm)) {
  if (!(arm %in% names(genes_by_arm))) {
    genes_by_arm[[arm]] <- character(0)
  }
  if (length(genes_by_arm[[arm]]) < 200) {
    remaining_genes <- commonGenes[!commonGenes %in% genes_by_arm[[arm]]]
    top_arm_genes <- remaining_genes[order(averageExpression[commonGenes %in% remaining_genes], decreasing = TRUE)[1:200]]
    genes_by_arm[[arm]] <- unique(c(genes_by_arm[[arm]], top_arm_genes))
  }
}

final_selected_genes <- unlist(genes_by_arm)
rawCounts <- rawCounts[final_selected_genes, , drop = FALSE]

normCounts <- log2(1 + rawCounts)
normCounts <- scale(normCounts, scale = FALSE)

if (scale_on_reference_label) {
  if (is.list(referenceCells)) {
    scaleFactor <- list()
    for (label in names(referenceCells)) {
      scaleFactor[[label]] <- rowMeans(normCounts[, referenceCells[[label]], drop = FALSE])
    }
    scaleFactor <- do.call(rbind, scaleFactor)
    scaleFactor <- apply(scaleFactor, 2, median)
  } else {
    scaleFactor <- rowMeans(normCounts[, referenceCells, drop = FALSE])
  }
} else {
  scaleFactor <- rowMeans(normCounts)
}

normCounts <- normCounts - scaleFactor
normCounts <- funTrim(normCounts, lo = -3, up = 3)

geneMetadata2 <- geneMetadata2[which(geneMetadata2$hgnc_symbol %in% topExprGenes), ]
geneMetadata2 <- geneMetadata2[order(geneMetadata2$chromosome_num, geneMetadata2$start_position), ]

genomicWindows <- lapply(c(1:23), function(chrom) {
  genesC <- geneMetadata2[which(geneMetadata2$chromosome_num == chrom), ]
  chr_arms <- unique(genesC$chr_arm)
  chrom_windows <- list()
  for (arm in chr_arms) {
    genesArm <- genesC[which(genesC$chr_arm == arm), ]
    N <- nrow(genesArm)
    iter <- round(window_size / 2)
    if (N > window_size) {
      gw <- lapply(seq(iter + 1, N - iter, by = window_step), function(i) {
        genesArm[(i - iter):(i + iter), "hgnc_symbol"] |> unlist() |> as.character()
      })
      names(gw) <- paste0(chrom, ".", arm, 1:length(gw))
    } else if (N > 0) {
      gw <- list(as.character(unlist(genesArm$hgnc_symbol)))
      names(gw) <- paste0(chrom, ".", arm, 1)
    } else {
      gw <- list()
    }
    chrom_windows <- c(chrom_windows, gw)
  }
  chrom_windows
})
genomicWindows <- unlist(genomicWindows, recursive = FALSE)

genomicScores <- sapply(genomicWindows, function(g) {
  if (length(g) == 1) {
    normCounts[g, ]
  } else {
    colMeans(normCounts[g, , drop = FALSE])
  }
})

if (scale_on_reference_label) {
  if (is.list(referenceCells)) {
    genomicScoresReferenceLabel <- do.call(rbind, lapply(referenceCells, function(cells) genomicScores[cells, , drop = FALSE]))
  } else {
    genomicScoresReferenceLabel <- genomicScores[referenceCells, , drop = FALSE]
  }
  Q01Q99 <- apply(genomicScoresReferenceLabel, 2, stats::quantile, probs = c(0 + threshold_percentile, 1 - threshold_percentile))
  genomicScoresTrimmed <- apply(genomicScores, 1, function(v) {
    v[which(v >= Q01Q99[1, ] & v <= Q01Q99[2, ])] <- 0
    v
  })
} else {
  if (is.na(reference_var) || is.null(obs)) {
    Q01Q99 <- apply(genomicScores, 2, stats::quantile, probs = c(0 + threshold_percentile, 1 - threshold_percentile))
    genomicScoresTrimmed <- apply(genomicScores, 1, function(v) {
      v[which(v >= Q01Q99[1, ] & v <= Q01Q99[2, ])] <- 0
      v
    })
  } else {
    cellLines <- split(colnames(rawCounts), obs[colnames(rawCounts), reference_var])
    high_threshold <- median(unlist(sapply(cellLines, function(z) apply(genomicScores[z, , drop = FALSE], 2, function(v) quantile(v, probs = c(1 - threshold_percentile))))))
    low_threshold <- median(unlist(sapply(cellLines, function(z) apply(genomicScores[z, , drop = FALSE], 2, function(v) quantile(v, probs = c(0 + threshold_percentile))))))
    genomicScoresTrimmed <- t(apply(genomicScores, 2, function(z) {
      z[which(z > low_threshold & z < high_threshold)] <- 0
      z
    }))
  }
}

raw_out <- t(as.matrix(genomicScores))
trim_out <- as.matrix(genomicScoresTrimmed)

write.table(raw_out, file.path(output_dir, "rawGenomicScores.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(trim_out, file.path(output_dir, "genomicScores.tsv"), sep = "\t", quote = FALSE, col.names = NA)
cnv_fraction <- colMeans(abs(trim_out) > 0)
write.table(data.frame(cnv_fraction = cnv_fraction), file.path(output_dir, "cell_metadata.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(data.frame(window = names(genomicWindows), n_genes = lengths(genomicWindows)), file.path(output_dir, "genomic_windows.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

