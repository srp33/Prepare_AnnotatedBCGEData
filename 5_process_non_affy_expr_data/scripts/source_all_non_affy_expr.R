library(GEOquery)
library(tidyverse)
library(janitor)
options(download.file.method.GEOquery = "wget")
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 1000) #helps with the download for GSE62944 which is a very large file

# create download directory for temporary files
tmp_dir <- "tmp/"
dir.create(tmp_dir)

# create data directory for saving data
data_dir <- "/Data/expression_data/"

source("scripts/GSE81538_expr.R")
source("scripts/GSE96058_expr.R")
source("scripts/GSE62944_expr.R")
source("scripts/ICGC_South_Korea_expr.R")
source("scripts/Metabric_expr.R")
source("scripts/SCAN_B_expr.R")

# delete temporary download directory
unlink(tmp_dir, recursive = TRUE, force = TRUE)