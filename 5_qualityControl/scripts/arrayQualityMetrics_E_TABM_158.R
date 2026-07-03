base_dir <- "E-TABM-158"
unlink(base_dir, recursive = TRUE)
if (!dir.exists(base_dir)) {
  dir.create(base_dir)
}

out_file_path <- "/Data/arrayQualityMetrics_results/E_TABM_158.tsv"

if (file.exists(out_file_path)) {
  print(paste0(out_file_path, " already exists."))
} else {
  download.file("https://zenodo.org/records/20097812/files/E-TABM-158.zip?download=1",
                destfile = paste0(base_dir, "/E_TABM_158.zip"), method = "wget")

  unzip(paste0(base_dir, "/E_TABM_158.zip"), exdir = base_dir)

  celFilePaths <- list.files(base_dir, pattern = "*.CEL", full.names = T, ignore.case = T)
print(celFilePaths)

  cel_files <- ReadAffy(filenames = celFilePaths)
print(cel_files)

  test_results <- arrayQualityMetrics(expressionset = cel_files, force = TRUE, outdir = str_c(base_dir, "/arrayQualityMetrics"))

  test_results$arrayTable %>%
    mutate(sampleNames = extractSampleIDs(sampleNames)) %>%
    mutate(heatmap_outlier = test_results$modules$heatmap@outliers@statistic > test_results$modules$heatmap@outliers@threshold) %>%
    mutate(boxplot_outlier = test_results$modules$boxplot@outliers@statistic > test_results$modules$boxplot@outliers@threshold) %>%
    mutate(maplot_outlier = test_results$modules$maplot@outliers@statistic > test_results$modules$maplot@outliers@threshold) %>%
    dplyr::select(sampleNames, heatmap_outlier, boxplot_outlier, maplot_outlier, any_of("ScanDate")) %>%
    write_tsv(out_file_path)

  print(paste("Saved to ", out_file_path))

  unlink(base_dir, recursive = TRUE, force = TRUE)
}
