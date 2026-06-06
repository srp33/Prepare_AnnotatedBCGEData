cancerTypeSamples <- getGEOSuppFiles("GSE62944", makeDirectory = F, baseDir = tmp_dir, filter_regex = "GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt.gz")
CancerType <- rownames(cancerTypeSamples) %>%
  read_tsv(col_names = F)
colnames(CancerType) <- c("Sample_ID", "cancer_type")
cancer_sample_ids <- CancerType %>%
  dplyr::filter(cancer_type == "BRCA") %>%
  dplyr::mutate(Sample_ID = substr(Sample_ID, 1, 12)) %>%
  pull(Sample_ID)

normalTypeSamples <- getGEOSuppFiles("GSE62944", makeDirectory = F, baseDir = tmp_dir, filter_regex = "GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt.gz")
normalType <- rownames(normalTypeSamples) %>%
  read_tsv(col_names = F)
colnames(normalType) <- c("Sample_ID", "cancer_type")
normal_sample_ids <- normalType %>%
  dplyr::filter(cancer_type == "BRCA") %>%
  mutate(Sample_ID = substr(Sample_ID, 1, 12)) %>%
  pull(Sample_ID)

# create a vector to replace the unknown variables with NA
na_strings <- c("NA", "[Unknown]", "[Not Available]", "[Not Evaluated]", "[Not Applicable]")
clinicalVariables <- getGEOSuppFiles("GSE62944", makeDirectory = F, baseDir = tmp_dir, filter_regex = "GSE62944_06_01_15_TCGA_24_548_Clinical_Variables_9264_Samples.txt.gz")
Clinical_Variables <- rownames(clinicalVariables) %>%
  read_tsv(col_names = F, na = na_strings) %>%
  dplyr::select(-(X2:X3))
Clinical_Variables[1,1] <- "specimen_id"
Clinical_Variables <- t(Clinical_Variables)
Clinical_Variables <- as_tibble(row_to_names(Clinical_Variables, 1, remove_row = TRUE)) %>%
  mutate(bcr_patient_barcode = substr(specimen_id, 1, 12)) %>%
  dplyr::rename(Sample_ID = bcr_patient_barcode) %>%
  dplyr::select(-specimen_id) %>%
  dplyr::filter(!is.na(bcr_patient_uuid)) %>%
  mutate(Dataset_ID = "GSE62944_Tumor") %>% 
  mutate(Platform_ID = "GPL9052") %>%
  dplyr::select(Dataset_ID, Sample_ID, Platform_ID, everything())

cancer_clinical <- filter(Clinical_Variables, Sample_ID %in% cancer_sample_ids) %>%
  dplyr::select(where(~ mean(is.na(.)) <= 0.5)) %>%
  dplyr::select(-form_completion_date, -prospective_collection, -retrospective_collection, -lymph_nodes_examined_count, -days_to_initial_pathologic_diagnosis, -icd_o_3_site, -informed_consent_verified, -patient_id, -birth_days_to, -tumor_tissue_site)
normal_clinical <- filter(Clinical_Variables, Sample_ID %in% normal_sample_ids) %>%
  dplyr::select(where(~ mean(is.na(.)) <= 0.5)) %>%
  dplyr::select(-form_completion_date, -prospective_collection, -retrospective_collection, -lymph_nodes_examined_count, -days_to_initial_pathologic_diagnosis, -icd_o_3_site, -informed_consent_verified, -patient_id, -birth_days_to, -tumor_tissue_site)

write_tsv(cancer_clinical, paste0(raw_metadata_dir, "GSE62944_Tumor.tsv"))
write_tsv(normal_clinical, paste0(raw_metadata_dir, "GSE62944_Normal.tsv"))

# summarise metadata variables
varSummary_tumor <- summariseVariables(cancer_clinical)
varSummary_normal <- summariseVariables(normal_clinical)

write_tsv(varSummary_tumor$numSummary, file.path("/Data/metadata_summaries/GSE62944_Tumor_num.tsv"))
write_tsv(varSummary_tumor$charSummary, file.path("/Data/metadata_summaries/GSE62944_Tumor_char.tsv"))

write_tsv(varSummary_normal$numSummary, file.path("/Data/metadata_summaries/GSE62944_Normal_num.tsv"))
write_tsv(varSummary_normal$charSummary, file.path("/Data/metadata_summaries/GSE62944_Normal_char.tsv"))

print("Writing GSE62944 to file")
write_tsv(cancer_clinical, paste0(data_dir, "GSE62944_Tumor.tsv"))
write_tsv(normal_clinical, paste0(data_dir, "GSE62944_Normal.tsv"))
