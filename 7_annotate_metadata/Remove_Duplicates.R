library(tidyverse)
library(readxl)
library(writexl)

############################################################################################
# This part of the code is not fully reproducible because it involves an interative process
# of checking for unmapped or inconsistently mapped metadata, revising the mappings and
# then checking again.
############################################################################################

categorical_initial = read_xlsx("Metadata_Mappings_Initial.xlsx", sheet = "Categorical", col_types = rep("text", 6)) %>%
  dplyr::mutate(primitive_type = "Categorical")

numeric_initial = read_xlsx("Metadata_Mappings_Initial.xlsx", sheet = "Numeric", col_types = rep("text", 6)) %>%
  dplyr::mutate(primitive_type = "Numeric")

categorical_initial = categorical_initial[!duplicated(categorical_initial), ]
numeric_initial = numeric_initial[!duplicated(numeric_initial), ]

write_xlsx(list("Categorical" = categorical_initial, "Numeric" = numeric_initial), "Metadata_Mappings_Deduplicated.xlsx")