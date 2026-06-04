# Identify directory where input variable (normalized data) is stored
datadir <- "/Data/expression_data"
metadata_dir <- "/Data/prelim_metadata"

# Turn a dataset into an ExpressionSet
makeExprs <- function(file_path, metadata) {
  GSEdata <- read_tsv(file_path) %>%
    rename_with(~str_replace_all(., "_.+", ""))

  gene_names <- pull(GSEdata, Gene)
  GSEdata <- dplyr::select(GSEdata, -Gene)
  data_matrix <- as.matrix(GSEdata)
  rownames(data_matrix) <- gene_names

  if (is.null(metadata)) {
    data_expr <- ExpressionSet(assayData = data_matrix)
  } else {
    sampleIDs <- dplyr::pull(metadata, Sample_ID)
    metadata <- dplyr::select(metadata, -Sample_ID) %>%
      as.data.frame()
    rownames(metadata) <- sampleIDs
    metadata <- AnnotatedDataFrame(metadata)

    data_expr <- ExpressionSet(assayData = data_matrix, phenoData = metadata)
  }
}

# Run dopplegangR using pairwise comparisons of datasets
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
        result <- doppelgangR(dopple_data, automatic.smokingguns = TRUE, BPPARAM = MulticoreParam(workers = 16))

        result_summary <- summary(result)
        write_tsv(result_summary, out_file_path)
    }

    done_this_time <- c(done_this_time, out_file_path)
  }
}

unlink("cache", recursive = TRUE, force = TRUE)
