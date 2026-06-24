# FYI: The first parts of this code are what we used when the files were available
#      on the ICGC portal. We are keeping it here for provenance. We placed a copy
#      of the processed files online and are using those.

# download.file("https://dcc.icgc.org/api/v1/download?fn=/current/Projects/BRCA-KR/donor.BRCA-KR.tsv.gz",
#               destfile = paste0(tmp_dir, "donor.BRCA-KR.tsv.gz"))
# donor_BRCA_KR <- read_tsv(paste0(tmp_dir, "donor.BRCA-KR.tsv.gz")) %>%
#   dplyr::select(-c(project_code, study_donor_involved_in, submitted_donor_id, donor_tumour_stage_at_diagnosis_supplemental, prior_malignancy,
#                    cancer_type_prior_malignancy, cancer_history_first_degree_relative))
# 
# 
# download.file("https://dcc.icgc.org/api/v1/download?fn=/current/Projects/BRCA-KR/specimen.BRCA-KR.tsv.gz",
#               destfile = paste0(tmp_dir, "specimen.BRCA-KR.tsv.gz"))
# specimen_BRCA_KR <- read_tsv(paste0(tmp_dir, "specimen.BRCA-KR.tsv.gz")) %>%
#   dplyr::filter(specimen_type == "Primary tumour - solid tissue") %>%
#   dplyr::select(-c(icgc_specimen_id, project_code, study_specimen_involved_in, submitted_specimen_id, submitted_donor_id, specimen_type_other,
#                    specimen_interval, specimen_donor_treatment_type_other, specimen_processing, specimen_processing_other, specimen_storage,
#                    specimen_storage_other, tumour_confirmed, specimen_biobank_id, specimen_available, tumour_grade_supplemental,
#                    tumour_stage_supplemental, digital_image_of_stained_section, percentage_cellularity, specimen_type, 
#                    specimen_donor_treatment_type, specimen_biobank))

# combined_df <- donor_BRCA_KR %>%
#   inner_join(specimen_BRCA_KR, by = "icgc_donor_id") %>%
#   dplyr::rename(Sample_ID = icgc_donor_id) %>% 
#   dplyr::select(-c("donor_diagnosis_icd10", "donor_tumour_staging_system_at_diagnosis")) %>%
#   mutate(Dataset_ID = "ICGC_KR", .before = Sample_ID) %>%
#   mutate(Platform_ID = "Illumina HiSeq", .after = Sample_ID)

tmp_file_path = paste0(tmp_dir, "ICGC.tsv.gz")

download.file("https://zenodo.org/records/20097812/files/ICGC_KR_sample_metadata.tsv.gz?download=1", destfile = tmp_file_path)

combined_df <- read_tsv(tmp_file_path) %>%
  mutate(Dataset_ID = "ICGC_KR", .before = Sample_ID) %>%
  mutate(Platform_ID = "GPL16791", .after = Sample_ID) %>%
  dplyr::rename(submitted_donor_id = submitted_donor_id.x) %>%
  dplyr::select(-project_code.x,
                -study_donor_involved_in,
                -donor_diagnosis_icd10,
                -donor_tumour_stage_at_diagnosis_supplemental,
                -donor_tumour_staging_system_at_diagnosis,
                -prior_malignancy,
                -cancer_type_prior_malignancy,
                -cancer_history_first_degree_relative,
                -icgc_specimen_id,
                -project_code.y,
                -study_specimen_involved_in,
                -submitted_specimen_id,
                -submitted_donor_id.y,
                -specimen_type,
                -specimen_type_other,
                -specimen_interval,
                -specimen_donor_treatment_type,
                -specimen_donor_treatment_type_other,
                -specimen_processing,
                -specimen_processing_other,
                -specimen_storage,
                -specimen_storage_other,
                -tumour_confirmed,
                -specimen_biobank,
                -specimen_biobank_id,
                -specimen_available,
                -tumour_grade_supplemental,
                -tumour_stage_supplemental,
                -digital_image_of_stained_section,
                -percentage_cellularity) %>%
  distinct()

#summarise metadata variables
varSummary <- summariseVariables(combined_df)

if (nrow(varSummary$numSummary) >= 1) {
  write_tsv(varSummary$numSummary, file.path("/Data/metadata_summaries/ICGC_KR_num.tsv"))
}

if (nrow(varSummary$charSummary) >= 1) {
  write_tsv(varSummary$charSummary, file.path("/Data/metadata_summaries/ICGC_KR_char.tsv"))
}

print("Writing ICGC_KR to file!")
write_tsv(combined_df, file.path(raw_metadata_dir, "ICGC_KR.tsv"))
write_tsv(combined_df, paste0(data_dir, "ICGC_KR.tsv"))

# d_BRCA_KR <- read_tsv(paste0(tmp_dir, "donor.BRCA-KR.tsv.gz")) 
# s_BRCA_KR <- read_tsv(paste0(tmp_dir, "specimen.BRCA-KR.tsv.gz"))

# combined_data <- d_BRCA_KR %>%
#   inner_join(s_BRCA_KR, by = "icgc_donor_id") %>%
#   dplyr::rename(Sample_ID = icgc_donor_id)
