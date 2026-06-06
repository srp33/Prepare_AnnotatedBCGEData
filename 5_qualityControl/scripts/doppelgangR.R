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

processCombo <- function(file_path1, file_path2, dataset_id1, dataset_id2, metadata_file_path1, metadata_file_path2, out_file_path) {
  if (file.exists(out_file_path)) {
    print(paste0(out_file_path, " already exists"))
  } else {
    metadata1 <- read_tsv(metadata_file_path1) %>%
      dplyr::select(-Dataset_ID, -Platform_ID)
    metadata2 <- read_tsv(metadata_file_path2) %>%
      dplyr::select(-Dataset_ID, -Platform_ID)

    if (ncol(metadata1) == 0 | ncol(metadata2) == 0) {
      metadata1 <- NULL
      metadata2 <- NULL
    }

    print(paste0("Reading in ", file_path1, "!"))
    print(paste0("Reading in ", file_path2, "!"))

    eSet1 <- makeExprs(file_path1, metadata1)
    eSet2 <- makeExprs(file_path2, metadata2)

    dopple_data <- list(eSet1, eSet2)
    names(dopple_data) <- c(dataset_id1, dataset_id2)

    result <- tryCatch(
    {
      message("Attempting with default for phenoFinder.args.")
      doppelgangR(
        dopple_data,
        automatic.smokingguns = TRUE,
        BPPARAM = SerialParam()#MulticoreParam(workers = 16)
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
        BPPARAM = SerialParam()#MulticoreParam(workers = 16)
      )
    })

    write_tsv(summary(result), out_file_path)
  }
}

# Enable parallelization
num_parallel = 8
registerDoParallel(num_parallel)
stopifnot(foreach::getDoParWorkers() == num_parallel)

# Run dopplegangR for pairwise comparisons of datasets.
file_paths <- list.files(datadir, full.names = T)

# Get all unique pairs.
pairs <- as.data.frame(
  t(combn(file_paths, 2)),
  stringsAsFactors = FALSE
)

set.seed(0)

colnames(pairs) <- c("file_path1", "file_path2")
pairs <- mutate(pairs, dataset_id1 = basename(file_path1)) %>%
  mutate(dataset_id2 = basename(file_path2)) %>%
  mutate(dataset_id1 = str_replace(dataset_id1, "\\.tsv\\.gz", "")) %>%
  mutate(dataset_id2 = str_replace(dataset_id2, "\\.tsv\\.gz", "")) %>%
  mutate(metadata_file_path1 = str_c(metadata_dir, "/", dataset_id1, ".tsv")) %>%
  mutate(metadata_file_path2 = str_c(metadata_dir, "/", dataset_id2, ".tsv")) %>%
  mutate(out_file_path = str_c("/Data/doppelgangR_results/", dataset_id1, "_", dataset_id2, ".tsv.gz")) %>%
  slice_sample(prop = 1)

foreach (i = 1:nrow(pairs)) %dopar% {
  row <- as.vector(as.matrix(pairs[i,]))
  file_path1 <- row[1]
  file_path2 <- row[2]
  dataset_id1 <- row[3]
  dataset_id2 <- row[4]
  metadata_file_path1 <- row[5]
  metadata_file_path2 <- row[6]
  out_file_path <- row[7]

  processCombo(file_path1, file_path2, dataset_id1, dataset_id2, metadata_file_path1, metadata_file_path2, out_file_path)
}

unlink("cache", recursive = TRUE, force = TRUE)
