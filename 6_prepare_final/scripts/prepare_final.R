library(biomaRt)
library(tidyverse)
library(tools)

# set biomart cache location. Without this, code fails becase of docker permissions issues.
Sys.setenv(BIOMART_CACHE = "/Data/cache")

normalized_data <- "/Data/expression_data"

IQRray_result <- "/Data/IQRray_results"

clean_colnames_expr_data <- "/Data/clean_colnames_expr_data/"
if (!dir.exists(clean_colnames_expr_data)) {
  dir.create(clean_colnames_expr_data)
}

IQRray_filtered <- "/Data/IQRray_filtered_data/"
if (!dir.exists(IQRray_filtered)) {
  dir.create(IQRray_filtered)
}

final_meta <- "/Data/analysis_ready_metadata/"
if (!dir.exists(final_meta)) {
  dir.create(final_meta)
}

meta_expr_matched_data <- "/Data/matched_expr_data/"
if (!dir.exists(meta_expr_matched_data)) {
  dir.create(meta_expr_matched_data)
}

analysis_ready_data <- "/Data/analysis_ready_expression_data/"
if (!dir.exists(analysis_ready_data)) {
  dir.create(analysis_ready_data)
}

source("scripts/clean_expression_data_colnames.R")
source("scripts/IQRray_filter_samples.R")
source("scripts/match_data.R")
source("scripts/add_gene_symbol.R")
source("scripts/merge_metadata_summaries.R")

unlink(clean_colnames_expr_data, recursive = TRUE, force = TRUE)
unlink(IQRray_filtered, recursive = TRUE, force = TRUE)
unlink(meta_expr_matched_data, recursive = TRUE, force = TRUE)
unlink(meta_dir, recursive = TRUE, force = TRUE)
unlink("/Data/cache", recursive = TRUE, force = TRUE) 
