
huExon <- read_tsv(paste0(IQRray_result, "/huExon.tsv")) %>%
  mutate(Platform = "HuEx") %>%
  mutate(Passing = TRUE)

huGene <- read_tsv(paste0(IQRray_result, "/huGene.tsv")) %>%
  mutate(Platform = "HuGene") %>%
  mutate(Passing = TRUE)

U95_2 <- read_tsv(paste0(IQRray_result, "/U95_2.tsv")) %>%
  mutate(Platform = "U95v2") %>%
  mutate(Passing = value > 45759.6835)

U133_A_Early_Access <- read_tsv(paste0(IQRray_result, "/U133A_Early_Access.tsv")) %>%
  mutate(Platform = "U133A") %>%
  mutate(Passing = value > 53812.1375)

U133_A <- read_tsv(paste0(IQRray_result, "/U133_A.tsv")) %>%
  mutate(Platform = "U133A") %>%
  mutate(Passing = value > 53812.1375)

U133_A_multiple_chip <- read_tsv(paste0(IQRray_result, "/U133_A_multiple_chip.tsv")) %>%
  mutate(Platform = "U133A") %>%
  mutate(Passing = value > 53812.1375)

U133_A2 <- read_tsv(paste0(IQRray_result, "/U133_A2.tsv")) %>%
  mutate(Platform = "U133A2") %>%
  mutate(Passing = value > 66795.84023)

U133_B <- read_tsv(paste0(IQRray_result, "/U133_B.tsv")) %>%
  mutate(Platform = "U133B") %>%
  mutate(Passing = value > 50700.16364)

U133_plus_2 <- read_tsv(paste0(IQRray_result, "/U133_plus_2.tsv")) %>%
  mutate(Platform = "U133Plus2") %>%
  mutate(Passing = value > 136266.0795)

U133_plus_2_multiple_chip <- read_tsv(paste0(IQRray_result, "/U133_Plus2_multiple_chip.tsv")) %>%
  mutate(Platform = "U133Plus2") %>%
  mutate(Passing = value > 136266.0795)

E_TABM_158 <- read_tsv(paste0(IQRray_result, "/E_TABM_158.tsv")) %>%
  mutate(Platform = "U133A") %>%
  mutate(Passing = value > 53812.1375)

big_list <- do.call("rbind", list(huExon, huGene, U95_2, U133_A_Early_Access, U133_A, U133_A_multiple_chip, U133_A2, U133_B, U133_plus_2, U133_plus_2_multiple_chip, E_TABM_158))

goodQuality <- big_list[which(big_list$Passing == "TRUE"), ]
goodQuality$gsmID = gsub("gsm", "GSM", goodQuality$gsmID)

goodQual <- function(expr_file) {
    cat("\n")
    print(paste0("Reading in ", expr_file, "!"))
    cat("\n")

    data <- read_tsv(expr_file)
    Sample_ID <- colnames(data)
    all_samples <- goodQuality %>%
      pull(gsmID)
    keep_samples <- intersect(Sample_ID, all_samples)
    clean_data <- dplyr::select(data, Dataset, Gene, all_of(keep_samples))
}

special_cases <- c("GSE62944_Tumor.tsv", "GSE62944_Normal.tsv", "GSE81538.tsv", "GSE96058_HiSeq.tsv",
"GSE96058_NextSeq.tsv", "ICGC_KR.tsv", "METABRIC.tsv")

for (file in list.files(clean_colnames_expr_data, full.names = T)) {

  file_name <- file %>% basename() %>% file_path_sans_ext()
  out_file_path <- paste0(IQRray_filtered, file_name, ".gz")

  if (file_name %in% special_cases) {
    new_data <- read_tsv(file)
  } else {
    new_data <- goodQual(file)
  }

  print(paste0("Writing ", out_file_path, " to file"))
  
  write_tsv(new_data, out_file_path)
}