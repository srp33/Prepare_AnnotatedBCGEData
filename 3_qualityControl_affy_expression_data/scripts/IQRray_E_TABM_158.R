# IQRray for E_TABM_158

base_dir <- "E-TABM-158/"
unlink(base_dir, recursive = TRUE)
if (!dir.exists(base_dir)) {
  dir.create(base_dir)
}

download.file("https://osf.io/download/mg8kj",
              destfile = paste0(base_dir, "/E_TABM_158.zip"), method = "wget")

unzip(paste0(base_dir, "E_TABM_158.zip"), exdir = base_dir)

celFilePaths <- list.files(base_dir, pattern = "*.CEL", full.names = T, ignore.case = T)

my_data <- ReadAffy(filenames = celFilePaths)
IQR_score <- IQRray_affy(my_data) %>%
  as_tibble(rownames = "celfileID")

# code below matches sample names from metadata file to IQR_score
# download metadata file
download.file("https://www.ebi.ac.uk/arrayexpress/files/E-TABM-158/E-TABM-158.sdrf.txt",
              destfile = paste0(base_dir, "ETABM_158_meta.txt"))

ETABM_meta <- read_tsv(paste0(base_dir, "ETABM_158_meta.txt"))

meta_col <- ETABM_meta %>%
  dplyr::select(c("Array Data File", "Source Name")) %>%
  dplyr::rename(celfileID = `Array Data File`) %>%
  dplyr::rename(gsmID = `Source Name`)

joint_cols <- full_join(meta_col, IQR_score) %>%
  mutate(gseID = "E_TABM_158") %>%
  dplyr::select(gseID, gsmID, celfileID, value)


write_tsv(joint_cols, paste0(IQRray_file_path, "E_TABM_158.tsv"))
print("Saved to E_TABM_158.tsv")

unlink(base_dir, recursive = TRUE, force = TRUE)