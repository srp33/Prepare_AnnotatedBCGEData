
# identify directory where input variable (normalized data) is stored
datadir <- "/Data/expression_data"

# turn each dataset into an Expressionset
makeExprs <- function(file_path) {
  GSEdata <- read_tsv(file_path) %>%
    rename_with(~str_replace_all(., "_.+", ""))
  gene_names <- pull(GSEdata, Gene)
  GSEdata <- dplyr::select(GSEdata, -Gene)
  data_matrix <- as.matrix(GSEdata)
  rownames(data_matrix) <- gene_names
  data_expr <- ExpressionSet(assayData = data_matrix)
}

# run dopplegangR using pairwise comparisons of datasets
file_paths <- list.files(datadir, full.names = T)
for (i in 1:(length(file_paths) - 1)) {
  file_path1 <- file_paths[i]
  for (j in (i + 1):length(file_paths)) {
    file_path2 <- file_paths[j]

    file_names <- basename(c(file_path1, file_path2)) %>%
      file_path_sans_ext() %>%
      file_path_sans_ext()

     out_file_path <- paste0(doppel_dir,  paste0(file_names, collapse = "_"), ".tsv")

    if (file.exists(out_file_path)) {
      print(paste0(out_file_path, " already exists"))
    } else {
        cat("\n")
        print(paste0("Reading in ", file_path1, "!"))
        print(paste0("Reading in ", file_path2, "!"))
        cat("\n")

        exprs_file1 <- makeExprs(file_path1)
        exprs_file2 <- makeExprs(file_path2)
        dopple_file <- list(exprs_file1, exprs_file2)
        names(dopple_file) <- file_names

        # To make doppelgangR faster, users may run in parallel using the code snipet withh the "BPPARAM = MulticoreParam()" arguement.
        # Note that errors may occur on computer systems not equiped to handle parallelization.
        # Users should assign workers based on number of available cores
        # On computer systems not equiped for parallelization users should use the code snipet withh the "BPPARAM = SerialParam()" arguement

        # result <- doppelgangR(dopple_file, phenoFinder.args = NULL, BPPARAM = SerialParam())
        result <- doppelgangR(dopple_file, phenoFinder.args = NULL, BPPARAM = MulticoreParam(workers = 16))

        result_summary <- summary(result)
        write_tsv(result_summary, out_file_path)
    }
  }
}

#delete cache dir
unlink("cache", recursive = TRUE, force = TRUE)
