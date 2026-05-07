#An update to readr caused GEOquery to fail.
#This temporary workaround fixed it so leaving it here incase it breaks again in the future

#Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 1000)
#readr::local_edition(1)

library(GEOquery)
library(tidyverse)
library(stringi)
library(janitor)
library(rlist)

# Define the file path to the unprocessed metadata directory
raw_metadata_dir <- "/Data/raw_metadata/"

# Create the unprocessed metadata folder if it doesn't exist
if (!dir.exists(raw_metadata_dir)) {
  dir.create(raw_metadata_dir, recursive = TRUE)
}

# Define the file path to the metadata directory
metadata_dir <- "/Data/prelim_metadata/"

# Create the metadata folder if it doesn't exist
if (!dir.exists(metadata_dir)) {
  dir.create(metadata_dir, recursive = TRUE)
}

# Define the file path to the variable summaries directory in tsv format
metadata_summaries <- "/Data/metadata_summaries/"

# Create the folder if it doesn't exist
if (!dir.exists(metadata_summaries)) {
  dir.create(metadata_summaries, recursive = TRUE)
}

#getting required functions
source("functions/getFromGEO.R")
source("functions/removeUnusefulCols.R")
source("functions/summariseVariables.R")
source("functions/writeOutput.R")

#parse Files
source("scripts/E_TABM_158_meta.R")
source("scripts/single_chips.R")
source("scripts/multiple_chips.R")
