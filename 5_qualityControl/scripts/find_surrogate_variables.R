datadir <- "/Data/expression_data4"
metadata_dir <- "/Data/prelim_metadata2"
out_dir <- "/Data/sva_results"

getExprData <- function(expr_file_path) {
  print(paste0("Reading ", expr_file_path))

  expr_data <- read_tsv(expr_file_path)
  gene_names <- pull(expr_data, Entrez_Gene_ID)
  expr_data <- dplyr::select(expr_data, -Dataset_ID, -Entrez_Gene_ID, -Ensembl_Gene_ID, -HGNC_Symbol, -Chromosomal_Band) %>%
    as.matrix()
  rownames(expr_data) <- gene_names

  return(expr_data)
}

getMetadata <- function(metadata_file_path) {
  print(paste0("Reading ", metadata_file_path))

  metadata <- read_tsv(metadata_file_path)# %>%
#    dplyr::select(-Dataset_ID)

  #metadata <- mutate(metadata, Platform_ID = str_c(Dataset_ID, "____", Platform_ID))
  platform_id <- pull(metadata, Platform_ID) %>%
    unique()

  metadata <- dplyr::select(metadata, -Platform_ID)
metadata <- dplyr::select(metadata, -Dataset_ID)

  if (ncol(metadata) == 0) {
    return(list(Platform_ID = platform_id, Metadata = NULL))
  } else {
    sampleIDs <- dplyr::pull(metadata, Sample_ID)
    metadata <- dplyr::select(metadata, -Sample_ID) %>%
      as.data.frame()
    rownames(metadata) <- sampleIDs

    return(list(Platform_ID = platform_id, Metadata = metadata))
  }
}

processDataset <- function(expr_file_path, metadata_file_path, is_microarray, out_file_path) {
  if (file.exists(out_file_path)) {
    return(NULL)
  }

#  expr_data <- getExprData(expr_file_path)
  metadata <- getMetadata(metadata_file_path)
#print(basename(metadata_file_path))
return(metadata$Platform_ID)
return(NULL)
#stop()
#  is_microarray <- dataset_id == "METABRIC" | dataset_id %in% affy_gseIDs
#
#  print(expr_file_path)
#  print(is_microarray)

  sample_ids <- sort(intersect(colnames(expr_data), rownames(metadata)))

  if (length(sample_ids) < 5) {
    stop(paste0("There are few, if any, matching samples between the metadata and expression data for ", file_path1, "."))
  }

  expr_data <- expr_data[,sample_ids]

  if (is.null(metadata)) {
    metadata <- metadata[sample_ids, , drop = FALSE]
    eSet <- ExpressionSet(assayData = expr_data, phenoData = AnnotatedDataFrame(metadata))
  } else {
    eSet <- ExpressionSet(assayData = expr_data)
  }

    #write_tsv(result, out_file_path)
}

# Enable parallelization
#num_parallel = 16
#registerDoParallel(num_parallel)
#stopifnot(foreach::getDoParWorkers() == num_parallel)

# Run dopplegangR for pairwise comparisons of datasets.
expr_file_paths <- list.files(datadir, full.names = T)

tmps <- c()
#foreach (i = 1:nrow(expr_file_paths)) %dopar% {
for (i in 1:length(expr_file_paths)) {
  dataset_id <- sub(".tsv.gz", "", basename(expr_file_paths[i]))

  expr_file_path <- expr_file_paths[i]
  metadata_file_path <- str_c(metadata_dir, "/", dataset_id, ".tsv")
  out_file_path <- str_c(out_dir, "/", dataset_id, ".tsv.gz")

  tmp <- processDataset(expr_file_path, metadata_file_path, is_microarray, out_file_path)
  tmps <- c(tmps, tmp)
#break
}

print(sort(unique(tmps)))
