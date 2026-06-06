# Fetch RAW data from GEO
if (!file.exists(paste0(tmp_dir, "GSE62944_RAW.tar"))) {
  GSE <- getGEOSuppFiles(GEO = "GSE62944", makeDirectory = F, baseDir = tmp_dir, filter_regex = "GSE62944_RAW.tar")
  storage_dir <- rownames(GSE)
  untar(storage_dir[1], exdir = tmp_dir)
}

tpm_file_path <- paste0(tmp_dir, "GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_TPM.txt.gz")
print(paste0("Reading ", tpm_file_path))
GSE62944_tumor_df <- read_tsv(tpm_file_path, col_names = F)

#Process tumor gene expression data
print("Reading tumor samples file")
cancerTypeSamples <- getGEOSuppFiles("GSE62944", makeDirectory = F, baseDir = tmp_dir, filter_regex = "GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt.gz")
CancerType <- rownames(cancerTypeSamples) %>%
  read_tsv(col_names = F)
colnames(CancerType) <- c("Sample_ID", "cancer_type")
CancerType <- CancerType %>%
  dplyr::filter(cancer_type == "BRCA")

print("Transposing tumor data")
Transposed_tumor_expr <- as.data.frame(t(GSE62944_tumor_df), stringsAsFactors = F)
Transposed_tumor_expr[1, 1] <- "Sample_ID"
Transposed_tumor_expr <- row_to_names(Transposed_tumor_expr, 1, remove_row = TRUE, remove_rows_above = TRUE)

print("Merging tumor data")
Merged_tumor_df <- Transposed_tumor_expr %>%
  inner_join(CancerType, by = "Sample_ID") %>%
  dplyr::filter(cancer_type == "BRCA") %>%
  dplyr::select(-cancer_type)

print("Transposing merged data")
GSE62944_tumor_data <- t(Merged_tumor_df) %>%
  row_to_names(1, remove_row = TRUE, remove_rows_above = TRUE) %>%
  as_tibble(rownames = "HGNC_Symbol")

sampleIDs <- names(GSE62944_tumor_data[2:ncol(GSE62944_tumor_data)])
sampleIDs <- substr(sampleIDs, 1, 12)
colnames(GSE62944_tumor_data)[2:ncol(GSE62944_tumor_data)] <- sampleIDs

# match with metadata file for tumor samples
print("Reading tumor metadata file")
tumor_metadata <- read_tsv("/Data/prelim_metadata/GSE62944_Tumor.tsv")

print("Creating filtered tibble")
mapping_tbl_tumor <- tibble(Expr_ID = names(GSE62944_tumor_data)[2:ncol(GSE62944_tumor_data)]) %>%
  filter(Expr_ID %in% tumor_metadata$Sample_ID)

print("Selecting columns")
clean_tumor_data <- dplyr::select(GSE62944_tumor_data, HGNC_Symbol, all_of(mapping_tbl_tumor$Expr_ID))


# Process normal gene expression data
print("Reading normal samples file")
normalTypeSamples <- getGEOSuppFiles("GSE62944", makeDirectory = F, baseDir = tmp_dir, filter_regex = "GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt.gz")
normalType <- rownames(normalTypeSamples) %>%
  read_tsv(col_names = F)
colnames(normalType) <- c("Sample_ID", "cancer_type")
normalType <- normalType %>%
  dplyr::filter(cancer_type == "BRCA")

print("Reading normal expression data")
GSE62944_normal_df <- read_tsv(paste0(tmp_dir, "GSM1697009_06_01_15_TCGA_24.normal_Rsubread_TPM.txt.gz"), col_names = F)

print("Transposing normal data")
Transposed_normal <- as.data.frame(t(GSE62944_normal_df), stringsAsFactors = F)
Transposed_normal[1, 1] <- "Sample_ID"
Transposed_normal <- row_to_names(Transposed_normal, 1, remove_row = TRUE, remove_rows_above = TRUE)

print("Merging normal data")
Merged_df_normal <- Transposed_normal %>%
  inner_join(normalType, by = "Sample_ID") %>%
  dplyr::filter(cancer_type == "BRCA") %>%
  dplyr::select(-cancer_type)

print("Transposing merged normal data")
GSE62944_normal_data <- t(Merged_df_normal) %>%
  row_to_names(1, remove_row = TRUE, remove_rows_above = TRUE) %>%
  as_tibble(rownames = "HGNC_Symbol")

sampleIDs <- names(GSE62944_normal_data[2:ncol(GSE62944_normal_data)])
sampleIDs <- substr(sampleIDs, 1, 12)
colnames(GSE62944_normal_data)[2:ncol(GSE62944_normal_data)] <- sampleIDs

# match with metadata file for normal samples
print("Reading normal metadata")
normal_metadata <- read_tsv("/Data/prelim_metadata/GSE62944_Normal.tsv")

print("Creating normal tibble")
mapping_tbl <- tibble(Expr_ID = names(GSE62944_normal_data)[2:ncol(GSE62944_normal_data)]) %>%
  filter(Expr_ID %in% normal_metadata$Sample_ID)

print("Selecting normal data columns")
clean_normal_data <- dplyr::select(GSE62944_normal_data, HGNC_Symbol, all_of(mapping_tbl$Expr_ID))

out_tumor_file_path <- paste0(data_dir, "GSE62944_Tumor.tsv.gz")
out_normal_file_path <- paste0(data_dir, "GSE62944_Normal.tsv.gz")

print(paste0("Writing GSE62944 to ", out_tumor_file_path, " and ", out_normal_file_path, "!"))
write_tsv(clean_tumor_data, out_tumor_file_path)
write_tsv(clean_normal_data, out_normal_file_path)
