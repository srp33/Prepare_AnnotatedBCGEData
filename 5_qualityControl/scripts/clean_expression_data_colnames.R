file_paths <- list.files("/Data/expression_data", full.names = T)

hgnc_special_cases <- c("GSE62944_Tumor", "GSE62944_Normal", "GSE81538", "GSE96058_HiSeq",
                        "GSE96058_NextSeq", "ICGC_KR", "METABRIC")

ensg_special_cases <- c("ABiM.100", "ABiM.405", "Normal.66", "OSLO2EMIT0.103", "SCANB.9206")

# remove extra characters in column names
for (file_path in file_paths) {
  out_file_path <- paste0("/Data/expression_data2/", basename(file_path))

  if (file.exists(out_file_path)) {
    print(paste0(out_file_path, " has already been processed."))
    next
  }

  cat("\n")
  print(paste0("Reading in ", file_path, "!"))
  cat("\n")

  gseID <- basename(file_path)
  gseID <- gsub(".tsv.gz", "", gseID)

  if (gseID %in% hgnc_special_cases) {
    expr_data <- read_tsv(file_path) %>%
      mutate(Dataset_ID = gseID, .before = everything())
  } else {
    if (gseID %in% ensg_special_cases) {
      expr_data <- read_tsv(file_path) %>%
        dplyr::rename(ENSG = Gene) %>%
        mutate(Dataset_ID = gseID, .before = everything())
    } else {
      expr_data <- read_tsv(file_path)

      names(expr_data) <- gsub("_.+", "", names(expr_data))
      names(expr_data) <- gsub("gsm", "GSM", names(expr_data))

      expr_data <- mutate(expr_data, Dataset_ID = gseID, .before = everything())
    }
  }

  write_tsv(expr_data, out_file_path)
}
