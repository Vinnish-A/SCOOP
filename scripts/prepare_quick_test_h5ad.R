#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(qs)
  library(scCustomize)
  library(Seurat)
  library(SeuratObject)
})

set.seed(20260614)

write_counts <- function(seu, dest, assay = "RNA") {
  scCustomize::as.anndata(
    x = seu,
    file_path = dirname(dest),
    file_name = basename(dest),
    assay = assay,
    main_layer = "counts",
    other_layers = NULL,
    transer_dimreduc = TRUE
  )
}

postprocess_h5ad <- function(dest) {
  py_path <- sprintf("r'''%s'''", normalizePath(dest, mustWork = TRUE))
  script <- sprintf(paste(
    "import anndata as ad",
    "adata = ad.read_h5ad(%s)",
    "adata.layers['counts'] = adata.X.copy()",
    "for col in ['nCount_RNA', 'nFeature_RNA']:",
    "    if col in adata.obs:",
    "        del adata.obs[col]",
    "adata.uns['schema_version'] = 'scsp_agent_sop.h5ad.v1'",
    "adata.write_h5ad(%s)",
    sep = "\n"
  ), py_path, py_path)
  tf <- tempfile(fileext = ".py")
  writeLines(script, tf)
  on.exit(unlink(tf), add = TRUE)
  status <- system2(Sys.which("python"), tf)
  if (!identical(status, 0L)) {
    stop("Failed to postprocess H5AD: ", dest)
  }
}

pick_cells_by_sample <- function(meta, sample_col, n_samples = 10, cells_per_sample = 1000) {
  sample_counts <- sort(table(meta[[sample_col]]), decreasing = TRUE)
  eligible <- sample_counts[sample_counts >= cells_per_sample]
  if (length(eligible) < n_samples) {
    stop(
      "Only ", length(eligible), " samples have at least ", cells_per_sample,
      " cells in column ", sample_col, "."
    )
  }

  selected_samples <- names(utils::head(eligible, n_samples))
  selected_cells <- unlist(
    lapply(selected_samples, function(sample_id) {
      cells <- rownames(meta)[meta[[sample_col]] == sample_id]
      sample(cells, cells_per_sample)
    }),
    use.names = FALSE
  )

  list(cells = selected_cells, samples = selected_samples, counts = eligible[selected_samples])
}

pick_cells_balanced_by_group <- function(meta, sample_col, group_col,
                                         samples_per_group = 8,
                                         cells_per_sample = 3000) {
  group_values <- if (identical(group_col, ".sample_prefix")) {
    sub("_.*", "", as.character(meta[[sample_col]]))
  } else {
    as.character(meta[[group_col]])
  }

  sample_group_counts <- aggregate(
    rep(1, nrow(meta)),
    by = list(sample_id = meta[[sample_col]], group_label = group_values),
    FUN = sum
  )
  colnames(sample_group_counts)[3] <- "source_cells"
  sample_group_counts <- sample_group_counts[order(
    sample_group_counts$group_label,
    -sample_group_counts$source_cells
  ), ]

  selected <- do.call(rbind, lapply(sort(unique(sample_group_counts$group_label)), function(group_id) {
    eligible <- sample_group_counts[
      sample_group_counts$group_label == group_id &
        sample_group_counts$source_cells >= cells_per_sample,
    ]
    if (nrow(eligible) < samples_per_group) {
      stop(
        "Only ", nrow(eligible), " samples in ", group_id, " have at least ",
        cells_per_sample, " cells."
      )
    }
    utils::head(eligible, samples_per_group)
  }))

  selected_cells <- unlist(
    lapply(selected$sample_id, function(sample_id) {
      cells <- rownames(meta)[meta[[sample_col]] == sample_id]
      sample(cells, cells_per_sample)
    }),
    use.names = FALSE
  )

  list(cells = selected_cells, selected = selected)
}

minimal_metadata <- function(seu, sample_col, dataset_label, source_dataset_col = NULL) {
  md <- seu@meta.data

  output <- data.frame(
    sample_id = as.character(md[[sample_col]]),
    batch_id = as.character(md[[sample_col]]),
    assay = "RNA",
    dataset = dataset_label,
    total_counts = as.numeric(md[["nCount_RNA"]]),
    n_genes_by_counts = as.numeric(md[["nFeature_RNA"]]),
    stringsAsFactors = FALSE,
    row.names = rownames(md)
  )

  if (!is.null(source_dataset_col)) {
    if (identical(source_dataset_col, ".sample_prefix")) {
      output$sample_group <- sub("_.*", "", as.character(md[[sample_col]]))
    } else if (source_dataset_col %in% colnames(md)) {
      output$source_dataset <- as.character(md[[source_dataset_col]])
    }
  }

  if ("percent.mt" %in% colnames(md)) {
    output$pct_counts_mt <- as.numeric(md[["percent.mt"]])
  } else if ("pct_mito" %in% colnames(md)) {
    output$pct_counts_mt <- as.numeric(md[["pct_mito"]])
  }

  if ("percent.ribo" %in% colnames(md)) {
    output$pct_counts_ribo <- as.numeric(md[["percent.ribo"]])
  } else if ("percent.RPS" %in% colnames(md)) {
    output$pct_counts_ribo <- as.numeric(md[["percent.RPS"]])
  } else if ("pct_ribo" %in% colnames(md)) {
    output$pct_counts_ribo <- as.numeric(md[["pct_ribo"]])
  }

  if ("percent.hb" %in% colnames(md)) {
    output$pct_counts_hb <- as.numeric(md[["percent.hb"]])
  }

  seu@meta.data <- output
  seu
}

drop_empty_features <- function(seu, assay = "RNA") {
  counts <- GetAssayData(seu, assay = assay, layer = "counts")
  keep_features <- rownames(counts)[Matrix::rowSums(counts) > 0]
  subset(seu, features = keep_features)
}

prepare_one <- function(source, dest, sample_col, dataset_label,
                        n_samples = 10, cells_per_sample = 1000) {
  message("Reading: ", source)
  seu <- qread(source, nthreads = 4)
  DefaultAssay(seu) <- "RNA"

  picked <- pick_cells_by_sample(
    meta = seu@meta.data,
    sample_col = sample_col,
    n_samples = n_samples,
    cells_per_sample = cells_per_sample
  )

  message("Selected samples:")
  print(data.frame(
    sample_id = picked$samples,
    source_cells = as.integer(picked$counts),
    sampled_cells = cells_per_sample,
    row.names = NULL
  ))

  seu <- subset(seu, cells = picked$cells)
  seu <- minimal_metadata(seu, sample_col = sample_col, dataset_label = dataset_label)
  seu <- drop_empty_features(seu)

  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest)) {
    file.remove(dest)
  }

  message("Writing: ", dest)
  write_counts(seu, dest)
  postprocess_h5ad(dest)
  message("Done: ", dest)
}

prepare_one_balanced_public <- function(source, dest, sample_col, group_col,
                                        dataset_label,
                                        samples_per_group = 8,
                                        cells_per_sample = 3000,
                                        balance_label = "Selected group balance") {
  message("Reading: ", source)
  seu <- qread(source, nthreads = 4)
  DefaultAssay(seu) <- "RNA"

  picked <- pick_cells_balanced_by_group(
    meta = seu@meta.data,
    sample_col = sample_col,
    group_col = group_col,
    samples_per_group = samples_per_group,
    cells_per_sample = cells_per_sample
  )

  message("Selected samples:")
  print(transform(picked$selected, sampled_cells = cells_per_sample), row.names = FALSE)
  message(balance_label, ":")
  print(aggregate(sample_id ~ group_label, data = picked$selected, FUN = length))

  seu <- subset(seu, cells = picked$cells)
  seu <- minimal_metadata(
    seu,
    sample_col = sample_col,
    dataset_label = dataset_label,
    source_dataset_col = group_col
  )
  seu <- drop_empty_features(seu)

  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest)) {
    file.remove(dest)
  }

  message("Writing: ", dest)
  write_counts(seu, dest)
  postprocess_h5ad(dest)
  message("Done: ", dest)
}

out_dir <- file.path(getwd(), "h5ad", "canonical", "quick_test")
target <- commandArgs(trailingOnly = TRUE)
target <- if (length(target) == 0) "all" else target[[1]]

require_env_path <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    stop("Set environment variable ", name, " to the source .qs file path.", call. = FALSE)
  }
  value
}

if (target %in% c("all", "public")) {
  prepare_one_balanced_public(
    source = require_env_path("SCOOP_PUBLIC_QS"),
    dest = file.path(out_dir, "public_O_GSE154795_24samples_3000cells_balanced.h5ad"),
    sample_col = "sample",
    group_col = "dataset",
    dataset_label = "public_O_GSE154795",
    samples_per_group = 8,
    cells_per_sample = 3000
  )
}

if (target == "all") {
  gc()
}

if (target %in% c("all", "internal")) {
  prepare_one_balanced_public(
    source = require_env_path("SCOOP_INTERNAL_QS"),
    dest = file.path(out_dir, "internal_overall_sim_24samples_3000cells_balanced.h5ad"),
    sample_col = "orig.ident",
    group_col = ".sample_prefix",
    dataset_label = "internal_overall_sim",
    samples_per_group = 8,
    cells_per_sample = 3000,
    balance_label = "Selected sample group balance"
  )
}
