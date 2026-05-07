
#fetch RAW file from GEO and store in dataDir
GSE <- getGEOSuppFiles(GEO = "GSE62944", makeDirectory = F, baseDir = tmp_dir, filter_regex = "GSE62944_RAW.tar")

# #unzip the tar file for access to internal files
storage_dir <- rownames(GSE)
untar(storage_dir[1], exdir = tmp_dir)

GSE62944_tumor_df <- read_tsv(paste0(tmp_dir, "GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_TPM.txt.gz"), col_names = F)

#Process tumor gene expression data
cancerTypeSamples <- getGEOSuppFiles("GSE62944", makeDirectory = F, baseDir = tmp_dir, filter_regex = "GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt.gz")
CancerType <- rownames(cancerTypeSamples) %>%
  read_tsv(col_names = F)
colnames(CancerType) <- c("Sample_ID", "cancer_type")
CancerType <- CancerType %>%
  dplyr::filter(cancer_type == "BRCA")

#rearrange tumor dataframe
Transposed_tumor_expr <- as.data.frame(t(GSE62944_tumor_df), stringsAsFactors = F)
Transposed_tumor_expr[1, 1] <- "Sample_ID"
Transposed_tumor_expr <- row_to_names(Transposed_tumor_expr, 1, remove_row = TRUE, remove_rows_above = TRUE)

Merged_tumor_df <- Transposed_tumor_expr %>%
  inner_join(CancerType, by = "Sample_ID") %>%
  dplyr::filter(cancer_type == "BRCA") %>%
  dplyr::select(-cancer_type)

GSE62944_tumor_data <- t(Merged_tumor_df) %>%
  row_to_names(1, remove_row = TRUE, remove_rows_above = TRUE) %>%
  as_tibble(rownames = "HGNC_Symbol")

# match with metadata file for tumor samples
tumor_metadata <- read_tsv("/Data/prelim_metadata/GSE62944_Tumor.tsv")

mapping_tbl_tumor <- tibble(Expr_ID = names(GSE62944_tumor_data)[2:ncol(GSE62944_tumor_data)]) %>%
  filter(Expr_ID %in% tumor_metadata$Sample_ID)

clean_tumor_data <- dplyr::select(GSE62944_tumor_data, HGNC_Symbol, all_of(mapping_tbl_tumor$Expr_ID))


# Process normal gene expression data
normalTypeSamples <- getGEOSuppFiles("GSE62944", makeDirectory = F, baseDir = tmp_dir, filter_regex = "GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt.gz")
normalType <- rownames(normalTypeSamples) %>%
  read_tsv(col_names = F)
colnames(normalType) <- c("Sample_ID", "cancer_type")
normalType <- normalType %>%
  dplyr::filter(cancer_type == "BRCA")

GSE62944_normal_df <- read_tsv(paste0(tmp_dir, "GSM1697009_06_01_15_TCGA_24.normal_Rsubread_TPM.txt.gz"), col_names = F)

#rearrange normal dataframe
Transposed_normal <- as.data.frame(t(GSE62944_normal_df), stringsAsFactors = F)
Transposed_normal[1, 1] <- "Sample_ID"
Transposed_normal <- row_to_names(Transposed_normal, 1, remove_row = TRUE, remove_rows_above = TRUE)

Merged_df_normal <- Transposed_normal %>%
  inner_join(normalType, by = "Sample_ID") %>%
  dplyr::filter(cancer_type == "BRCA") %>%
  dplyr::select(-cancer_type)

GSE62944_normal_data <- t(Merged_df_normal) %>%
  row_to_names(1, remove_row = TRUE, remove_rows_above = TRUE) %>%
  as_tibble(rownames = "HGNC_Symbol")

# match with metadata file for normal samples
normal_metadata <- read_tsv("/Data/prelim_metadata/GSE62944_Normal.tsv")

mapping_tbl <- tibble(Expr_ID = names(GSE62944_normal_data)[2:ncol(GSE62944_normal_data)]) %>%
  mutate(Sample_ID = str_sub(Expr_ID, 1, 12)) %>%
  filter(Sample_ID %in% normal_metadata$bcr_patient_barcode)

clean_normal_data <- dplyr::select(GSE62944_normal_data, HGNC_Symbol, all_of(mapping_tbl$Expr_ID))


print("Writing GSE62944 to file!")
write_tsv(clean_tumor_data, paste0(data_dir, "GSE62944_Tumor.tsv.gz"))
write_tsv(clean_normal_data, paste0(data_dir, "GSE62944_Normal.tsv.gz"))
