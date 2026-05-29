library(tidyverse)
library(readxl)
library(writexl)

############################################################################################
# This part of the code is not fully reproducible because it involves an interative process
# of checking for unmapped or inconsistently mapped metadata, revising the mappings and
# then checking again.
############################################################################################

numeric_all = read_xlsx("Metadata_Mappings_Deduplicated.xlsx", sheet = "Numeric", col_types = rep("text", 7)) %>%
  separate_longer_delim(orig_values, delim = "||")

categorical_all = read_xlsx("Metadata_Mappings_Deduplicated.xlsx", sheet = "Categorical", col_types = rep("text", 7)) %>%
  separate_longer_delim(orig_values, delim = "||") %>%
  dplyr::select(dataset, orig_field, orig_values, NCIT_values)

categorical_temp = dplyr::select(categorical_all, dataset, orig_field)
numeric_temp = dplyr::select(numeric_all, dataset, orig_field)
common_temp = dplyr::intersect(categorical_temp, numeric_temp) %>%
  arrange(dataset, orig_field)