#clean up expression data names

file_paths <- list.files(normalized_data, full.names = T)

special_cases <- c("GSE62944_Tumor", "GSE62944_Normal", "GSE81538", "GSE96058_HiSeq",
                   "GSE96058_NextSeq", "ICGC_KR", "METABRIC")


# remove extra characters in column names
for (file in file_paths) {
  cat("\n")
  print(paste0("Reading in ", file, "!"))
  cat("\n")

  gseID <- basename(file)
  gseID <- gsub(".tsv.gz", "", gseID)

  if (gseID %in% special_cases) {
    expr_data <- read_tsv(file) %>%
      mutate(Dataset_ID = gseID, .before = everything())
  } else {
    expr_data <- read_tsv(file) %>%
      mutate(Dataset_ID = gseID, .before = everything()) %>%
      filter(!str_detect(Gene, "^AFFX"))
  }
  names(expr_data) <- gsub("_.+", "", names(expr_data))
  names(expr_data) <- gsub("gsm", "GSM", names(expr_data))
  
  write_tsv(expr_data, paste0(clean_colnames_expr_data, gseID, ".tsv.gz"))
}
