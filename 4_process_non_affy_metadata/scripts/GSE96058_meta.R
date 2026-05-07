
# get ExpressionSet from GEO for "this" GEO tag, create data frame of relevant information
gseID <- getGEO("GSE96058")

df_1 <- gseID[[1]]

metadata_1 <- pData(df_1)

# write un-curated metadata to file
write_tsv(metadata_1, file.path(raw_metadata_dir, "GSE96058_HiSeq.tsv"))

#pull phenotype data from GEO data frame
metadata_1 <- metadata_1 %>%
  clean_names() %>%
  removeCols() %>%
  rename_with(~str_replace_all(., "_ch1", "")) %>%
  dplyr::select(-c(title, description, scan_b_external_id, er_prediction_mgc, er_prediction_sgc, her2_prediction_mgc, her2_prediction_sgc,
                   ki67_prediction_mgc, ki67_prediction_sgc, nhg_prediction_mgc, pgr_prediction_mgc, pgr_prediction_sgc)) %>%
  dplyr::rename(Sample_ID = geo_accession) %>%
  mutate(Dataset_ID = "GSE96058_HiSeq", .before = Sample_ID) %>%
  mutate(Platform_ID = platform_id, .after = Sample_ID) %>%
  dplyr::select(-platform_id)

cols_to_change <- c("er_status", "pgr_status", "her2_status", "ki67_status")

# metadata_1 <- metadata_1 %>%
#   mutate_at(all_of(cols_to_change), ~ str_replace(., "0", "negative")) %>%
#   mutate_at(all_of(cols_to_change), ~ str_replace(., "1", "positive"))

#summarise metadata variables
varSummary <- summariseVariables(metadata_1)

if (nrow(varSummary$numSummary) >= 1) {
  write_tsv(varSummary$numSummary, file.path("/Data/metadata_summaries/GSE96058_HiSeq_num.tsv"))
}

if (nrow(varSummary$charSummary) >= 1) {
  write_tsv(varSummary$charSummary, file.path("/Data/metadata_summaries/GSE96058_HiSeq_char.tsv"))
}

df_2 <- gseID[[2]]

metadata_2 <- pData(df_2)

# write un-curated metadata to file
write_tsv(metadata_2, file.path(raw_metadata_dir, "GSE96058_NextSeq.tsv"))

#pull phenotype data from GEO data frame
metadata_2 <- metadata_2 %>%
  clean_names() %>%
  removeCols() %>%
  rename_with(~str_replace_all(., "_ch1", "")) %>%
  dplyr::select(-c(title, description, scan_b_external_id, er_prediction_mgc, er_prediction_sgc, her2_prediction_mgc, her2_prediction_sgc,
                   ki67_prediction_mgc, ki67_prediction_sgc, nhg_prediction_mgc, pgr_prediction_mgc, pgr_prediction_sgc)) %>%
  dplyr::rename(Sample_ID = geo_accession) %>%
  mutate(Dataset_ID = "GSE96058_NextSeq", .before = Sample_ID) %>%
  mutate(Platform_ID = platform_id, .after = Sample_ID) %>%
  dplyr::select(-platform_id)

# metadata_2 <- metadata_2 %>%
#   mutate_at(all_of(cols_to_change), ~ str_replace(., "0", "negative")) %>%
#   mutate_at(all_of(cols_to_change), ~ str_replace(., "1", "positive"))

#summarise metadata variables
varSummary <- summariseVariables(metadata_2)

if (nrow(varSummary$numSummary) >= 1) {
  write_tsv(varSummary$numSummary, file.path("/Data/metadata_summaries/GSE96058_NextSeq_num.tsv"))
}

if (nrow(varSummary$charSummary) >= 1) {
  write_tsv(varSummary$charSummary, file.path("/Data/metadata_summaries/GSE96058_NextSeq_char.tsv"))
}

print("Writing GSE96058 to file!")
write_tsv(metadata_1, paste0(data_dir, "GSE96058_HiSeq.tsv"))
write_tsv(metadata_2, paste0(data_dir, "GSE96058_NextSeq.tsv"))
