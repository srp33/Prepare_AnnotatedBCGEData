library(doParallel)
library(tidyverse)
library(arrayQualityMetrics)

extractSampleIDs <- function(sampleNames) {
  sampleNames <- gsub("_.+", "", sampleNames)
  sampleNames <- gsub("gsm", "GSM", sampleNames)
  return(sampleNames)
}

processDataset <- function(gseID, use_oligo_annotations) {
  out_file_path <- str_c("/Data/arrayQualityMetrics_results/", gseID, ".tsv.gz")

  if (file.exists(out_file_path)) {
    print(str_c(out_file_path, " already exists."))
  } else {
    print(str_c("Applying arrayQualityMetrics to ", gseID))

    tmp_dir_path = str_c(tempdir(), "/", gseID)

    if (dir.exists(tmp_dir_path)) {
      unlink(tmp_dir_path)
    }

    dir.create(tmp_dir_path)

    GSE <- getGEOSuppFiles(gseID, makeDirectory = F, baseDir = tmp_dir_path, filter_regex = "RAW.tar$")
    tmp <- rownames(GSE)
    untar(tmp[1], exdir = tmp_dir_path)

    cel_file_paths <- list.files(tmp_dir_path, pattern = "*.CEL", full.names = T, ignore.case = T)

    if (use_oligo_annotations) {
      cel_files = read.celfiles(cel_file_paths)
    } else {
      cel_files <- ReadAffy(filenames = cel_file_paths)
    }

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

datasets <- read_tsv("Affymetrix_Platforms.tsv") %>%
  mutate(use_oligo_annotations = geneChip %in% c("Affymetrix Human Exon 1.0 ST Array [transcript (gene) version]", "Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]")) %>%
  filter(gseID != "GSE7378") %>% # This dataset results in an error or an unknown reason (internal code error).
  filter(gseID != "GSE5460") # This resulted in a memory allocation error. It also indicated that GSM125119 and GSM125120 CEL files are corrupted.

# Enable parallelization
num_parallel = 8
registerDoParallel(num_parallel)
stopifnot(foreach::getDoParWorkers() == num_parallel)

#for (i in 1:nrow(datasets)) {
foreach (i = 1:nrow(datasets)) %dopar% {
  dataset <- datasets[i,]
  gseID <- pull(dataset, gseID)
  use_oligo_annotations <- pull(dataset, use_oligo_annotations)

  processDataset(gseID, use_oligo_annotations)
}
