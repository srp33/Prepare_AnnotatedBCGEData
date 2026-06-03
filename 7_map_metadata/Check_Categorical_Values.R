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

identifier_rows = dplyr::filter(categorical_all, str_detect(NCIT_values, "C25364")) %>%
  dplyr::select(dataset, orig_field, orig_values)

test_dataset <- NULL
# test_dataset <- "GSE20194"

for (metadata_file_path in list.files(path = "../Data/prelim_metadata", full.names = TRUE)) {
  this_dataset = basename(metadata_file_path)
  this_dataset = sub(".tsv", "", this_dataset)

  if (!is.null(test_dataset) && this_dataset != test_dataset)
    next
  
  numeric_dataset = dplyr::filter(numeric_all, dataset == this_dataset)

  categorical_dataset = categorical_all %>%
    dplyr::filter(dataset == this_dataset)

  metadata_dataset = read_tsv(metadata_file_path) %>%
    dplyr::select(-Dataset_ID, -Sample_ID, -Platform_ID)

  all_fields_dataset <- NULL

  for (field in colnames(metadata_dataset)) {
    values = pull(metadata_dataset, all_of(field)) %>%
      unique() %>%
      as.character()
    values[is.na(values)] <- "NA"

    fields = tibble(dataset = this_dataset, orig_field = field, orig_values = values)
    all_fields_dataset = bind_rows(all_fields_dataset, fields)
  }

  if (is.null(categorical_dataset) | is.null(all_fields_dataset)) {
    next
  }
  
  all_fields_dataset <- arrange(all_fields_dataset, orig_field, orig_values)

  diff1_dataset <- anti_join(all_fields_dataset, categorical_dataset) %>%
    anti_join(identifier_rows, by = join_by(dataset, orig_field)) %>%
    anti_join(numeric_dataset, by = join_by(dataset, orig_field))

  diff2_dataset <- anti_join(categorical_dataset, all_fields_dataset) %>%
    anti_join(identifier_rows, by = join_by(dataset, orig_field)) %>%
    anti_join(numeric_dataset, by = join_by(dataset, orig_field)) %>%
    arrange(orig_field, orig_values)

  if (nrow(diff1_dataset) > 0) {
    View(diff1_dataset)
    stop(paste0("There are missing annotations in diff1_dataset for ", this_dataset, "."))
  }
  
  if (nrow(diff2_dataset) > 0) {
    all_fields_dataset <- group_by(all_fields_dataset, dataset, orig_field) %>%
      summarize(data_values = str_flatten(sort(orig_values), collapse = "||"))
    diff2_dataset <- left_join(diff2_dataset, all_fields_dataset)
    View(diff2_dataset)
    stop(paste0("There are extra annotations in diff2_dataset for ", this_dataset, "."))
  }
  
  print(paste0("Annotations for ", this_dataset, " appear to be complete."))
}