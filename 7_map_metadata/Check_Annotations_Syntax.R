library(tidyverse)
library(readxl)
library(writexl)

############################################################################################
# This part of the code is not fully reproducible because it involves an interative process
# of checking for unmapped or inconsistently mapped metadata, revising the mappings and
# then checking again.
############################################################################################

categorical_fields_need_fixing = read_xlsx("Metadata_Mappings_Deduplicated.xlsx", sheet = "Categorical", col_types = rep("text", 7)) %>%
  separate_longer_delim(NCIT_field, delim = "||") %>%
  dplyr::filter(str_detect(NCIT_field, pattern = " \\(Code C\\d+\\)$", negate = TRUE))

categorical_values_need_fixing = read_xlsx("Metadata_Mappings_Deduplicated.xlsx", sheet = "Categorical", col_types = rep("text", 7)) %>%
  separate_longer_delim(NCIT_values, delim = "||") %>%
  dplyr::filter(str_detect(NCIT_values, pattern = " \\(Code C\\d+\\)", negate = TRUE))

numeric_fields_need_fixing = read_xlsx("Metadata_Mappings_Deduplicated.xlsx", sheet = "Numeric", col_types = rep("text", 7)) %>%
  separate_longer_delim(NCIT_field, delim = "||") %>%
  dplyr::filter(str_detect(NCIT_field, pattern = " \\(Code C\\d+\\)$", negate = TRUE))

numeric_values_need_fixing = read_xlsx("Metadata_Mappings_Deduplicated.xlsx", sheet = "Numeric", col_types = rep("text", 7)) %>%
  separate_longer_delim(NCIT_values, delim = "||") %>%
  dplyr::filter(str_detect(NCIT_values, pattern = " \\(Code C\\d+\\)$", negate = TRUE))