
# save supplementary files to Directory and read into separate data frames
cancerTypeSamples <- getGEOSuppFiles("GSE62944", makeDirectory = F, baseDir = tmp_dir, filter_regex = "GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt.gz")
CancerType <- rownames(cancerTypeSamples) %>%
  read_tsv(col_names = F)
colnames(CancerType) <- c("Sample_ID", "cancer_type")
CancerType <- CancerType %>%
  dplyr::filter(cancer_type == "BRCA")

normalTypeSamples <- getGEOSuppFiles("GSE62944", makeDirectory = F, baseDir = tmp_dir, filter_regex = "GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt.gz")
normalType <- rownames(normalTypeSamples) %>%
  read_tsv(col_names = F)
colnames(normalType) <- c("Sample_ID", "cancer_type")
normalType <- normalType %>%
  dplyr::filter(cancer_type == "BRCA") %>%
  mutate(bcr_patient_barcode = substr(Sample_ID, 1, 12))

# create a vector to replace the unknown variables with NA
na_strings <- c("NA", "[Unknown]", "[Not Available]", "[Not Evaluated]", "[Not Applicable]")
clinicalVariables <- getGEOSuppFiles("GSE62944", makeDirectory = F, baseDir = tmp_dir, filter_regex = "GSE62944_06_01_15_TCGA_24_548_Clinical_Variables_9264_Samples.txt.gz")
Clinical_Variables <- rownames(clinicalVariables) %>%
  read_tsv(col_names = F, na = na_strings) %>%
  dplyr::select(- (X2:X3))

# rearrange Clinical_Variables
Transposed_df <- as_tibble(t(Clinical_Variables), stringsAsFactors = F)
Transposed_df[1, 1] <- "Sample_ID"
Transposed_df <- row_to_names(Transposed_df, 1, remove_row = TRUE, remove_rows_above = TRUE)

#merge tumor and normal data frames by Sampleid
Merged_tumor_df <- Transposed_df %>%
  inner_join(CancerType, by = "Sample_ID") %>%
  dplyr::select(-cancer_type)

Merged_normal_df <- Transposed_df %>%
  inner_join(normalType, by = "bcr_patient_barcode") %>%
  rename(Sample_ID = Sample_ID.x) %>%
  dplyr::select(-c(Sample_ID.y, cancer_type))

# write un-curated metadata to file
write_tsv(Merged_tumor_df, paste0(raw_metadata_dir, "GSE62944_Tumor.tsv"))
write_tsv(Merged_normal_df, paste0(raw_metadata_dir, "GSE62944_Normal.tsv"))

# remove samples without an ID, and remove samples with only one value
Transposed_df <- Transposed_df %>%
  dplyr::filter(bcr_patient_uuid != "NA") %>%
  dplyr::select(-c("form_completion_date", "prospective_collection", "retrospective_collection", "tissue_source_site", 
                   "days_to_initial_pathologic_diagnosis", "icd_o_3_site", "lymph_nodes_examined_count", "tumor_tissue_site")) %>%
  remove_constant()


#calculate percentage of columns with missing data for easy visualisation
# NA_values <- as_tibble(colMeans(is.na(Merged_tumor_df)), rownames = "variable") %>%
#   pivot_wider(names_from = variable, values_from = value)
# new_df <- rbind(NA_values, Merged_tumor_df)

# keep columns with less than 50% NA
tumor_filtered <- Merged_tumor_df %>%
  dplyr::filter(bcr_patient_uuid != "NA") %>%
  dplyr::select(-c("form_completion_date", "prospective_collection", "retrospective_collection", "tissue_source_site", 
                   "days_to_initial_pathologic_diagnosis", "icd_o_3_site", "lymph_nodes_examined_count", "tumor_tissue_site")) %>%
  remove_constant()

tumor_filtered <- tumor_filtered[, colMeans(is.na(tumor_filtered)) < 0.5] %>%
  mutate(Dataset_ID = "GSE62944_Tumor", .before = Sample_ID) %>%
  mutate(Platform_ID = "GPL9052", .after = Sample_ID)


normal_filtered <- Merged_normal_df %>%
  dplyr::filter(bcr_patient_uuid != "NA") %>%
  dplyr::select(-c("form_completion_date", "prospective_collection", "retrospective_collection", "tissue_source_site", 
                   "days_to_initial_pathologic_diagnosis", "icd_o_3_site", "lymph_nodes_examined_count", "tumor_tissue_site")) %>%
  remove_constant()

normal_filtered <- normal_filtered[, colMeans(is.na(normal_filtered)) < 0.5] %>%
  mutate(Dataset_ID = "GSE62944_Normal", .before = Sample_ID) %>%
  mutate(Platform_ID = "GPL9052", .after = Sample_ID)


# code to figure out differences in columns between tumor data and normal data
# ***her2_ihc_score is absent in normal_filtered***
# a <- tibble(sort(colnames(tumor_filtered)))
# b <- tibble(sort(colnames(normal_filtered)))
# new_row = c(1)
# b = rbind(b,new_row)
# c <- cbind(a,b)

# summarise metadata variables
varSummary_tumor <- summariseVariables(tumor_filtered)
varSummary_normal <- summariseVariables(normal_filtered)

write_tsv(varSummary_tumor$numSummary, file.path("/Data/metadata_summaries/GSE62944_Tumor_num.tsv"))
write_tsv(varSummary_tumor$charSummary, file.path("/Data/metadata_summaries/GSE62944_Tumor_char.tsv"))

write_tsv(varSummary_normal$numSummary, file.path("/Data/metadata_summaries/GSE62944_Normal_num.tsv"))
write_tsv(varSummary_normal$charSummary, file.path("/Data/metadata_summaries/GSE62944_Normal_char.tsv"))


print("Writing GSE62944 to file!")
write_tsv(tumor_filtered, paste0(data_dir, "GSE62944_Tumor.tsv"))
write_tsv(normal_filtered, paste0(data_dir, "GSE62944_Normal.tsv"))
