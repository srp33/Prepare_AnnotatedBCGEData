# This script removes columns from the metadata.

removeUnusefulCols <- function(metadata) {

  #these columns are unuseful
  metadata <- metadata |>
    dplyr::select(-(starts_with(c("contact", "relation", "supplementary_file", "data_processing", "extract_protocol"))))

  # other unuseful columns which are common to multiple data sets
  remove_cols <- c("biomaterial_provider_ch1", "channel_count", "data_row_count", "filename", "growth_protocol_ch1", "hyb_protocol", 
  "id_ch1", "label_ch1", "label_protocol_ch1", "last_update_date", "molecule_ch1", "organism_ch1", "processor_id", "samplename_ch1", 
  "sample_name_ch1", "scan_protocol", "source_name_ch1", "status", "submission_date", "taxid_ch1", "treatment_protocol_ch1", "type")

  for (element in remove_cols) {
    if (element %in% names(metadata)) {
      metadata <- metadata |>
      dplyr::select(-all_of(element))
    }
  }

  column_names_to_remove <- c()

  for (i in seq_along(metadata)) {
           
    columnName <- colnames(metadata)[i]
    if (!str_detect(columnName, "^characteristics_")) {
        next
    }

    # The dataset GSE45255 has information in this column which is useful, so we are keeping it
    if (str_detect(columnName, "characteristics_ch1_9") & ("GSE45255" %in% metadata$Dataset_ID)) {
      next
    }  

    # This dataset has some inconsistencies which make it difficult for GEOquery to parse the file properly, thus it is done manually
    if (str_detect(columnName, "characteristics_ch1_4") & ("GSE22093" %in% metadata$Dataset_ID)) {
      next
    }    

    aCol <- pull(metadata, columnName)
    aCol <- aCol[!is.na(aCol)]
    proportion_with_colon <- sum(str_detect(aCol, ": ?")) / length(aCol)

    if (proportion_with_colon > 0.5) {
        column_names_to_remove <- c(column_names_to_remove, columnName)
    }    
  }

  if (length(column_names_to_remove) > 0) {
    metadata <- metadata |>
      dplyr::select(-all_of(column_names_to_remove))
  }

  return(metadata)
}
