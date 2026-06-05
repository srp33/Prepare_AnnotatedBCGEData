datadir <- "/Data/expression_data4"
metadata_dir <- "/Data/prelim_metadata2"

# Turn a dataset into an ExpressionSet
makeExprs <- function(file_path, metadata) {
  expr_data <- read_tsv(file_path)

  gene_names <- pull(expr_data, Entrez_Gene_ID)
  expr_data <- dplyr::select(expr_data, -Dataset_ID, -Entrez_Gene_ID, -Ensembl_Gene_ID, -HGNC_Symbol, -Chromosomal_Band) %>%
      as.matrix()
  rownames(expr_data) <- gene_names

  if (is.null(metadata)) {
    eSet <- ExpressionSet(assayData = expr_data)
  } else {
    sampleIDs <- dplyr::pull(metadata, Sample_ID)
    metadata <- dplyr::select(metadata, -Sample_ID) %>%
      as.data.frame()
    rownames(metadata) <- sampleIDs

    sampleIDs <- intersect(sampleIDs, colnames(expr_data))

    if (length(sampleIDs) < 5) {
      stop("There are few, if any, matching samples between the metadata and expression data.")
    }

    metadata <- metadata[sampleIDs, , drop = FALSE]
    metadata <- AnnotatedDataFrame(metadata)

    expr_data <- expr_data[,sampleIDs]

    eSet <- ExpressionSet(assayData = expr_data, phenoData = metadata)
  }

  return(eSet)
}

# Run dopplegangR for pairwise comparisons of datasets.
file_paths <- list.files(datadir, full.names = T)

done_this_time <- c()

for (i in 1:length(file_paths)) {
  file_path1 <- file_paths[i]

  for (j in 1:length(file_paths)) {
    file_path2 <- file_paths[j]

    if (file_path1 == file_path2) {
      next
    }

    dataset_ids <- basename(c(file_path1, file_path2)) %>%
      file_path_sans_ext() %>%
      file_path_sans_ext()

    # By sorting, we make sure we don't have to process the same ones
    # with their names in reverse order.
    combined_names <- paste0(sort(dataset_ids), collapse = "_")
    out_file_path <- paste0(doppel_dir, combined_names, ".tsv")

    if (file.exists(out_file_path)) {
      if (!(out_file_path %in% done_this_time)) {
        print(paste0(out_file_path, " already exists"))
      }
    } else {
        metadata_file_path1 <- paste0(metadata_dir, "/", dataset_ids[1], ".tsv")
        metadata_file_path2 <- paste0(metadata_dir, "/", dataset_ids[2], ".tsv")

        metadata1 <- read_tsv(metadata_file_path1) %>%
          dplyr::select(-Dataset_ID, -Platform_ID)
        metadata2 <- read_tsv(metadata_file_path2) %>%
          dplyr::select(-Dataset_ID, -Platform_ID)

        if (ncol(metadata1) == 0 | ncol(metadata2) == 0) {
          metadata1 <- NULL
          metadata2 <- NULL
        }

        cat("\n")
        print(paste0("Reading in ", file_path1, "!"))
        print(paste0("Reading in ", file_path2, "!"))
        cat("\n")

        eSet1 <- makeExprs(file_path1, metadata1)
        eSet2 <- makeExprs(file_path2, metadata2)

        dopple_data <- list(eSet1, eSet2)
        names(dopple_data) <- dataset_ids

        # result <- doppelgangR(dopple_data, phenoFinder.args = NULL, BPPARAM = SerialParam())

        result <- tryCatch(
        {
          message("Attempting with default for phenoFinder.args.")
          doppelgangR(
            dopple_data,
            automatic.smokingguns = TRUE,
            BPPARAM = MulticoreParam(workers = 16)
          )
        },
        error = function(e) {
          message("doppelgangR failed with phenoFinder.args.")
          message("Original error: ", conditionMessage(e))
          message("Retrying with phenoFinder.args = NULL")

          doppelgangR(
            dopple_data,
            phenoFinder.args = NULL,
            automatic.smokingguns = TRUE,
            BPPARAM = MulticoreParam(workers = 16)
          )
        })

        write_tsv(summary(result), out_file_path)
    }

    done_this_time <- c(done_this_time, out_file_path)
  }
}

unlink("cache", recursive = TRUE, force = TRUE)
