# required libraries
library(tidyverse)
library(SCAN.UPC)
library(doParallel)
library(janitor)
library(GEOquery)
library(tools)

# Enable parallelization
num_parallel = 2
registerDoParallel(num_parallel)
stopifnot(foreach::getDoParWorkers() == num_parallel)

# Download may time out, due to large file sizes, use this option to allow for longer download times
options(timeout = max(3000, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

# Create the normalized data folder if it doesn't exist
normalized_data <- "/Data/expression_data/"
if (!dir.exists(normalized_data)) {
    dir.create(normalized_data)
}

source("scripts/filter_chips.R")
#source("scripts/normalize_E_TABM_158.R")
#source("scripts/normalize_GSE23720.R")
source("scripts/normalize_single_chips.R")
# source("scripts/normalize_multiple_chips.R")
