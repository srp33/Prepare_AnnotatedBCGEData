library(GEOquery)
library(tidyverse)
library(janitor)
library(readxl)

# this setting helps with the download for GSE62944 which is a very large file
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10000)

#create download directory for temporary files
tmp_dir <- "tmp/"

# Define the file path to the unprocessed metadata directory
raw_metadata_dir <- "/Data/raw_metadata/"

# Define the file path to the metadata directory for saving data
data_dir <- "/Data/prelim_metadata/"

meta_summaries_dir <- "/Data/metadata_summaries/"

create_directory <- function(directory_path) {
  if (!dir.exists(directory_path)) {
    dir.create(directory_path, recursive = TRUE)
  }
}

unlink(tmp_dir, recursive = TRUE, force = TRUE)
create_directory(tmp_dir)
create_directory(raw_metadata_dir)
create_directory(data_dir)
create_directory(meta_summaries_dir)

source("functions/summariseVariables.R")
source("functions/removeCols.R")

#source("scripts/GSE81538_meta.R")
#source("scripts/GSE96058_meta.R")
#source("scripts/GSE62944_meta.R")
#source("scripts/METABRIC_meta.R")
#source("scripts/ICGC_South_Korea_meta.R")
source("scripts/SCAN_B_meta.R")

unlink(tmp_dir, recursive = TRUE, force = TRUE)
