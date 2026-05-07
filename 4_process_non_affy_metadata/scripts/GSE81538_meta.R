# get ExpressionSet from GEO for "this" GEO tag, create data frame of relevant information
gseID <- getGEO("GSE81538")
df <- gseID[[1]]

metadata <- pData(df)

# write un-curated metadata to file
write_tsv(metadata, file.path(raw_metadata_dir, "GSE81538.tsv"))

metadata <- metadata |>
  clean_names() |>
  removeCols() |>
  dplyr::select(-c("title", "description", "tissue_ch1")) |>
  rename_with(~str_replace_all(., "_ch1", "")) |>
  dplyr::rename(Sample_ID = geo_accession) |>
  mutate(Dataset_ID = "GSE81538", .before = Sample_ID) |>
  mutate(Platform_ID = platform_id, .after = Sample_ID) |>
  dplyr::select(-platform_id)

#summarise metadata variables
varSummary <- summariseVariables(metadata)

if (nrow(varSummary$numSummary) >= 1) {
  write_tsv(varSummary$numSummary, file.path("/Data/metadata_summaries/GSE81538_num.tsv"))
}

if (nrow(varSummary$charSummary) >= 1) {
  write_tsv(varSummary$charSummary, file.path("/Data/metadata_summaries/GSE81538_char.tsv"))
}

out_file_path <- paste0(data_dir, "GSE81538.tsv")
print(paste0("Writing GSE81538 to ", out_file_path))
write_tsv(metadata, out_file_path)