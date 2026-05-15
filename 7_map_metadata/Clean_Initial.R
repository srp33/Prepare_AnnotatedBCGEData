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

all_initial = bind_rows(categorical_initial, numeric_initial)

# Remove duplicated rows.
all_initial = all_initial[!duplicated(all_initial), ]

all_fields_values = NULL

for (metadata_file_path in list.files(path = "../Data/prelim_metadata", full.names = TRUE)) {
  dataset = basename(metadata_file_path)
  dataset = sub(".tsv", "", dataset)
  
  metadata_for_dataset = read_tsv(metadata_file_path)

  metadata = dplyr::select(metadata_for_dataset, -Dataset_ID, -Sample_ID, -Platform_ID)

  for (field in colnames(metadata)) {
    values = pull(metadata, all_of(field)) %>%
      unique() %>%
      as.character()

    fields_for_dataset = tibble(dataset = dataset, orig_field = field, orig_values = values)
    all_fields_values = bind_rows(all_fields_values, fields_for_dataset)
  }
}

# Identify which columns are in the metadata but have not been mapped.
# We adjust the metadata processing steps to ignore metadata fields
# that don't make sense to use (and leave comments explaining them).
# We ensure the rest are mapped to ontology term(s).

all_initial_fields = dplyr::select(all_initial, dataset, orig_field)
all_fields = dplyr::select(all_fields_values, dataset, orig_field)

anti_join(all_fields, all_initial_fields) %>%
  distinct() %>%
  write_xlsx("Metadata_Unmapped.xlsx")