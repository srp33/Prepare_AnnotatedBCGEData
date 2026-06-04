removeCols <- function(metadata) {

  metadata <- metadata %>%
    dplyr::select(- (starts_with(c("contact", "sample name", "relation", "supplementary", "data_processing", "library", "description_"))))

  #other unuseful columns
  remove_cols <- c("biomaterial_provider_ch1", "channel_count", "data_processing", "data_row_count", "extract_protocol_ch1", "extract_protocol_ch1_1",
                    "growth_protocol_ch1", "hyb_protocol", "instrument_model", "instrument_model_ch1", "label_ch1", "label_protocol_ch1", "last_update_date",
                   "molecule_ch1", "organism_ch1", "processor_id", "scan_protocol", "source_name_ch1", "status",
                   "submission_date", "taxid_ch1", "treatment_protocol_ch1", "type")

  for (element in remove_cols) {
    if (element %in% names(metadata)) {
      metadata <- metadata %>%
      dplyr::select(-all_of(element))
    }
  }

  column_names_to_remove <- c()

  for (i in seq_along(metadata)) {
    columnName <- colnames(metadata)[i]
    if (!str_detect(columnName, "^characteristics_")) {
      next
    }

    aCol <- pull(metadata, columnName)
    # aCol <- aCol[!is.na(aCol)]
    proportion_with_colon <- sum(str_detect(aCol, ": ?")) / length(aCol)

    if (proportion_with_colon > 0.5) {
      column_names_to_remove <- c(column_names_to_remove, columnName)
    }
  }

  if (length(column_names_to_remove) > 0) {
    metadata <- metadata %>%
      dplyr::select(-all_of(column_names_to_remove))
  }

  return(metadata)
}