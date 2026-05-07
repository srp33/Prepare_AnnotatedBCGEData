# Check https://data.mendeley.com/datasets/yzxtxn4nmd/3 for potential updates to the supplemental data.
# This code used the most recent version at the time of writing: version 2023-01-13
source_file_url <- "https://data.mendeley.com/public-files/datasets/yzxtxn4nmd/files/33d08b30-4685-4814-ab87-606a20c3092b/file_downloaded"
source_file_name <- "metadata_table.xlsx"

meta_file_path <- paste0(tmp_dir, source_file_name)

platform_id <- "Illumina HumanHT-12 V4.0 expression beadchip"

# Fill out columns you want to exclude
excluded_columns <- list(
    ABiM.100 = c("gex_type", "library_protocol"),
    ABiM.405 = c("gex_type", "sample", "library_protocol"),
    Normal.66 = c("gex_type", "fraction_duplication", "aligned_pairs", "library_protocol", "rna_qc_rin_rqs"),
    OSLO2EMIT0.103 = c("gex_type", "library_protocol"),
    SCANB.9206 = c("gex_type", "follow_up_cohort", "test_set", "training_set", "merged", "library", "rna", "sample", "case", "patient",
    "fraction_duplication", "aligned_pairs", "pm_reads", "pt_reads", "pf_reads", "read_string", "sequencer_serial", "pool_name",
    "library_protocol", "library_barcode", "library_batch_no", "rna_qc_rin_rqs", "rna_nd_conc", "qiacube_batch_no", "reference_year",
    "reference_source", "sampling_date", "days_to_lab", "minutes_to_rna_later")
)

# Download excel file (contains multiple sheets)
download.file(url = source_file_url, destfile = meta_file_path)#, mode = "wb") #use the wb mode if working on windows

clean_metadata <- function(sheet, path) {
    data_file <- path |>
        read_excel(sheet = sheet, col_types = "text", .name_repair = "universal") |>
        mutate(across(everything(), ~type.convert(.x, as.is = TRUE))) |>
        write_tsv(paste0(raw_metadata_dir, sheet, ".tsv")) |> # write un-curated metadata to file
        clean_names() |>
        rename(Sample_ID = gex_assay) |>
        mutate(Dataset_ID = sheet, .before = Sample_ID) |>
        mutate(Platform_ID = platform_id, .after = Sample_ID) |>
        dplyr::select(-excluded_columns[[sheet]]) |>
        write_tsv(paste0(data_dir, sheet, ".tsv")) |>
        summariseVariables()
    if (nrow(data_file$numSummary) >= 1) {      
        write_tsv(data_file$numSummary, file.path(paste0(meta_summaries_dir, sheet, "_num.tsv")))            
    }
    if (nrow(data_file$charSummary) >= 1) {
        write_tsv(data_file$charSummary, file.path(paste0(meta_summaries_dir, sheet, "_char.tsv")))
    }
}

meta_file_path %>%
    excel_sheets() %>%
    map(clean_metadata, path = meta_file_path)