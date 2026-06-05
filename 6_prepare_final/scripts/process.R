library(tidyverse)
library(tools)

source("scripts/IQRray_filter_samples.R")
source("scripts/merge_metadata_summaries.R")

unlink(clean_colnames_expr_data, recursive = TRUE, force = TRUE)
unlink(IQRray_filtered, recursive = TRUE, force = TRUE)
unlink(meta_expr_matched_data, recursive = TRUE, force = TRUE)
unlink(meta_dir, recursive = TRUE, force = TRUE)
