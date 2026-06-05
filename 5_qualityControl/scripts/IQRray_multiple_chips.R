out_file_path <- "/Data/IQRray_results/U133_A_multiple_chip.tsv"

if (!file.exists(out_file_path)) {
  # process GSE1456
  GSE1456 <- getGEO("GSE1456")
  GSE1456_U133A_geo_accession <- pData(GSE1456[[1]]) %>%
    mutate(geo_number = parse_number(as.character(geo_accession))) %>%
    dplyr::select(geo_accession, geo_number)

  GSE1456_U133B_geo_accession <- pData(GSE1456[[2]]) %>%
    mutate(geo_number = parse_number(as.character(geo_accession))) %>%
    dplyr::select(geo_accession, geo_number)

  # #process GSE3494
  GSE3494 <- getGEO("GSE3494")
  GSE3494_U133A_geo_accession <- pData(GSE3494[[1]]) %>%
    mutate(geo_number = parse_number(as.character(geo_accession))) %>%
    dplyr::select(geo_accession, geo_number)

  GSE3494_U133B_geo_accession <- pData(GSE3494[[2]]) %>%
    mutate(geo_number = parse_number(as.character(geo_accession))) %>%
    dplyr::select(geo_accession, geo_number)

  # #process GSE4922
  GSE4922 <- getGEO("GSE4922")
  GSE4922_U133A_geo_accession <- pData(GSE4922[[1]]) %>%
    mutate(geo_number = parse_number(as.character(geo_accession))) %>%
    dplyr::select(geo_accession, geo_number)

  GSE4922_U133B_geo_accession <- pData(GSE4922[[2]]) %>%
    mutate(geo_number = parse_number(as.character(geo_accession))) %>%
    dplyr::select(geo_accession, geo_number)

  #process GSE6532
  GSE6532 <- getGEO("GSE6532")
  GSE6532_U133A_geo_accession <- pData(GSE6532[[2]]) %>%
    mutate(geo_number = parse_number(as.character(geo_accession))) %>%
    dplyr::select(geo_accession, geo_number)

  GSE6532_U133B_geo_accession <- pData(GSE6532[[3]]) %>%
    mutate(geo_number = parse_number(as.character(geo_accession))) %>%
    dplyr::select(geo_accession, geo_number)

  GSE6532_U133_Plus2_geo_accession <- pData(GSE6532[[1]]) %>%
    mutate(geo_number = parse_number(as.character(geo_accession))) %>%
    dplyr::select(geo_accession, geo_number)

  final_score_GSE1456_U133A <- run_IQRray_multiple_chips("GSE1456", GSE1456_U133A_geo_accession)
  final_score_GSE3494_U133A <- run_IQRray_multiple_chips("GSE3494", GSE3494_U133A_geo_accession)
  final_score_GSE4922_U133A <- run_IQRray_multiple_chips("GSE4922", GSE4922_U133A_geo_accession)
  final_score_GSE6532_U133A <- run_IQRray_multiple_chips("GSE6532", GSE6532_U133A_geo_accession)

  big_IQR_file_U133A <- do.call("rbind", list(final_score_GSE1456_U133A, final_score_GSE3494_U133A, final_score_GSE4922_U133A, final_score_GSE6532_U133A))
  write_tsv(big_IQR_file_U133A, out_file_path)
  print(paste0("Saved to ", out_file_path))
}

out_file_path <- "/Data/IQRray_results/U133_B_multiple_chip.tsv"

if (!file.exists(out_file_path)) {
  final_score_GSE1456_U133B <- run_IQRray_multiple_chips("GSE1456", GSE1456_U133B_geo_accession)
  final_score_GSE3494_U133B <- run_IQRray_multiple_chips("GSE3494", GSE3494_U133B_geo_accession)
  final_score_GSE4922_U133B <- run_IQRray_multiple_chips("GSE4922", GSE4922_U133B_geo_accession)
  final_score_GSE6532_U133B <- run_IQRray_multiple_chips("GSE6532", GSE6532_U133B_geo_accession)

  big_IQR_file_U133B <- do.call("rbind", list(final_score_GSE1456_U133B, final_score_GSE3494_U133B, final_score_GSE4922_U133B, final_score_GSE6532_U133B))
  write_tsv(big_IQR_file_U133B, out_file_path)
  print(paste0("Saved to ", out_file_path))
}

out_file_path <- "/Data/IQRray_results/U133_Plus2_multiple_chip.tsv"

if (!file.exists(out_file_path)) {
  final_score_GSE6532_U133_Plus2 <- run_IQRray_multiple_chips("GSE6532", GSE6532_U133_Plus2_geo_accession)
  write_tsv(final_score_GSE6532_U133_Plus2, out_file_path)
  print(paste0("Saved to ", out_file_path))
}

# delete temporary download directory
unlink("GSE1456", recursive = TRUE, force = TRUE)
unlink("GSE3494", recursive = TRUE, force = TRUE)
unlink("GSE4922", recursive = TRUE, force = TRUE)
unlink("GSE6532", recursive = TRUE, force = TRUE)
