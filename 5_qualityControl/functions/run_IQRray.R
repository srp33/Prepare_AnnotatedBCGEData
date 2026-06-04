
# function to run the IQRray function across multiple datasets

# IQray on single chips
run_IQRray <- function(gseID) {

  # Define the file path to a temp directory for saving RAW data
  tmp_dir <- paste0("/tmp/", gseID)
  dir.create(tmp_dir)

  print(paste0("Downloading ", gseID, " for processing!"))

  GSE <- getGEOSuppFiles(gseID, makeDirectory = F, baseDir = tmp_dir, filter_regex = "RAW.tar$")
  tmp <- rownames(GSE)
  untar(tmp[1], exdir = tmp_dir)

  celFilePaths <- list.files(tmp_dir, pattern = "*.CEL", full.names = T, ignore.case = T)

  if (gseID %in% oligo_arrays) {
    my_data <- read.celfiles(filenames = celFilePaths)
    IQR_score <- IQRray_oligo(my_data) %>%
      as_tibble(rownames = "celfileID")
  } else {
    my_data <- ReadAffy(filenames = celFilePaths)
    IQR_score <- IQRray_affy(my_data) %>%
      as_tibble(rownames = "celfileID")
  }

  IQR_score <- mutate(IQR_score, gsmID = str_extract(IQR_score$celfileID, "[A-Za-z]+\\d+")) %>%
    relocate(gsmID)
  IQR_score <- cbind(gseID, IQR_score)

  unlink(tmp_dir, recursive = TRUE, force = TRUE)
  return(IQR_score)
}

# IQray on multiple chips
run_IQRray_multiple_chips <- function(gseID, geo_accession) {

  # Define the file path to a temp directory for saving RAW data
  tmp_dir <- paste0("/tmp/", gseID)
  dir.create(tmp_dir)

  print(paste0("Downloading ", gseID, " for processing!"))

  GSE <- getGEOSuppFiles(gseID, makeDirectory = F, baseDir = tmp_dir, filter_regex = "RAW.tar$")
  tmp <- rownames(GSE)
  untar(tmp[1], exdir = tmp_dir)

  celFiles <- list.files(tmp_dir, pattern = "*.CEL", full.names = T, ignore.case = T) %>%
    as_tibble()

  celFiles_new <- celFiles %>%
    mutate(gse_ID = as.numeric(str_extract(value, "\\d+")),
           geo_number = as.numeric(str_extract(value, "\\d+(?=\\.[A-Za-z])")))

  celFilePaths <- inner_join(geo_accession, celFiles_new) %>%
    dplyr::select(value) %>%
    pull()

  my_data <- ReadAffy(filenames = celFilePaths)
  IQR_score <- IQRray_affy(my_data) %>%
    as_tibble(rownames = "celfileID")

  IQR_score <- mutate(IQR_score, gsmID = str_extract(IQR_score$celfileID, "[A-Za-z]+\\d+")) %>%
    relocate(gsmID)
  IQR_score <- cbind(gseID, IQR_score)

  unlink(tmp_dir, recursive = TRUE, force = TRUE)
  return(IQR_score)
}
