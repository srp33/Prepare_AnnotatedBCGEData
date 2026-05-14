out_file_path <- paste0(metadata_dir, "E_TABM_158.tsv")

if (file.exists(out_file_path)) {
  print(paste0(out_file_path, " already exists, so it will not be processed."))
} else {
  tmp_dir <- "tmp/"
  unlink(tmp_dir, recursive = TRUE, force = TRUE)
  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir)
  }

  #download metadata file
  download.file("https://www.ebi.ac.uk/arrayexpress/files/E-TABM-158/E-TABM-158.sdrf.txt",
                destfile = paste0(tmp_dir, "ETABM_158_meta.txt"))

  #read the sample and data relationship file into a table
  etabm_158 <- read_tsv(paste0(tmp_dir, "ETABM_158_meta.txt"))

  # write un-curated metadata to file
  write_tsv(etabm_158, file.path(raw_metadata_dir, "E_TABM_158.tsv"))

  #filter out unuseful columns
  SDRF <- etabm_158 %>%
    dplyr::select(c("Source Name", starts_with("Characteristics"))) %>%
    dplyr::select(-c("Characteristics [BioSourceType]", "Characteristics [DiseaseState]", "Characteristics [OrganismPart]",
                   "Characteristics [Organism]")) %>%
    clean_names() %>%
    rename(Sample_ID = source_name) %>%
    rename_with(~str_replace_all(., "characteristics_", "")) %>%
    mutate(Dataset_ID = "E_TABM_158", .before = Sample_ID) %>%
    mutate(Platform_ID = "GPL4685", .after = Sample_ID)

  SDRF <- replace(SDRF, SDRF=='n/a', NA)

  varSummary <- summariseVariables(SDRF)

  if (nrow(varSummary$numSummary) >= 1) {
    write_tsv(varSummary$numSummary, file.path(metadata_summaries, "E_TABM_158_num.tsv"))
  }

  if (nrow(varSummary$charSummary) >= 1) {
    write_tsv(varSummary$charSummary, file.path(metadata_summaries, "E_TABM_158_char.tsv"))
  }

  print("Writing E_TABM_158 to file!")
  write_tsv(SDRF, out_file_path)

  unlink(tmp_dir, recursive = TRUE, force = TRUE)
}
