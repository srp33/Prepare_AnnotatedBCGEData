extractSampleIDs <- function(sampleNames) {
  sampleNames <- gsub("_.+", "", sampleNames)
  sampleNames <- gsub("gsm", "GSM", sampleNames)
  return(sampleNames)
}

processDataset <- function(gseID, supp_data_index, suffix) {
  tmp_dir_path <- str_c(tempdir(), "/", gseID)
  out_file_path <- str_c("/Data/arrayQualityMetrics_results/", gseID, suffix, ".tsv.gz")

  if (file.exists(out_file_path)) {
    print(str_c(out_file_path, " already exists"))
  } else {
    if (dir.exists(tmp_dir_path)) {
      unlink(tmp_dir_path, recursive = TRUE, force = TRUE)
    }
    dir.create(tmp_dir_path)

    # Retrieve accession IDs for this subset.
    supp_file_data <- getGEO(gseID)

    geo_accessions <- pData(supp_file_data[[supp_data_index]]) %>%
      mutate(geo_number = parse_number(as.character(geo_accession))) %>%
      pull(geo_accession)

    print(paste0("Downloading expression data from ", gseID, " for processing"))

    GSE <- getGEOSuppFiles(gseID, makeDirectory = F, baseDir = tmp_dir_path, filter_regex = "RAW.tar$")
    tmp <- rownames(GSE)
    untar(tmp[1], exdir = tmp_dir_path)

    cel_file_paths <- list.files(tmp_dir_path, pattern = "*.CEL", full.names = T, ignore.case = T)
    cel_file_paths <- cel_file_paths[grepl(paste(geo_accessions, collapse = "|"), cel_file_paths, ignore.case = T)]

    cel_files <- ReadAffy(filenames = cel_file_paths)

    test_results <- arrayQualityMetrics(expressionset = cel_files, force = TRUE, outdir = str_c(tmp_dir_path, "/arrayQualityMetrics"))

    test_results$arrayTable %>%
      mutate(sampleNames = extractSampleIDs(sampleNames)) %>%
      mutate(heatmap_outlier = test_results$modules$heatmap@outliers@statistic > test_results$modules$heatmap@outliers@threshold) %>%
      mutate(boxplot_outlier = test_results$modules$boxplot@outliers@statistic > test_results$modules$boxplot@outliers@threshold) %>%
      mutate(maplot_outlier = test_results$modules$maplot@outliers@statistic > test_results$modules$maplot@outliers@threshold) %>%
      dplyr::select(sampleNames, heatmap_outlier, boxplot_outlier, maplot_outlier, any_of("ScanDate")) %>%
      write_tsv(out_file_path)

    unlink(tmp_dir_path, recursive = TRUE)
  }
}

processDataset("GSE1456", 1, "_U133A")
processDataset("GSE1456", 2, "_U133B")
processDataset("GSE3494", 1, "_U133A")
processDataset("GSE3494", 2, "_U133B")
processDataset("GSE4922", 1, "_U133A")
processDataset("GSE4922", 2, "_U133B")
processDataset("GSE6532", 2, "_U133A")
processDataset("GSE6532", 3, "_U133B")
processDataset("GSE6532", 1, "_U133Plus2")
