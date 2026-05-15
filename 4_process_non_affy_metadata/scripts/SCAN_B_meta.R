# Check https://data.mendeley.com/datasets/yzxtxn4nmd/3 for potential updates to the supplemental data.
# This code used the most recent version at the time of writing: version 2023-01-13
source_file_url <- "https://data.mendeley.com/public-files/datasets/yzxtxn4nmd/files/33d08b30-4685-4814-ab87-606a20c3092b/file_downloaded"
source_file_name <- "metadata_table.xlsx"

meta_file_path <- paste0(tmp_dir, source_file_name)

platform_id <- "Illumina HumanHT-12 V4.0 expression beadchip"

# Fill out columns you want to exclude
excluded_columns <- list(
    ABiM.100 = c("gex_type", "library_protocol", "clin_group", "eval_group", "ssp_cc15", "ssp_pam50", "ssp_subtype", "ssp_ror_as_t0", "ssp_ror_num_tadj", "ssp_ror_risk_cat", "ssp_ror_binary_risk_cat", "ssp_etr", "ssp_etr_dc", "prosigna_ffpe_subtype", "prosigna_ffpe_ror", "prosigna_ffpe_ror_risk_cat", "prosigna_ffpe_ror_binary_risk_cat", "prosigna_ffpe_etr", "prosigna_ffpe_etr_dc", "prosigna_ft_subtype", "prosigna_ft_ror", "prosigna_ft_ror_risk_cat", "prosigna_ft_ror_binary_risk_cat", "prosigna_ft_etr", "prosigna_ft_etr_dc"),
    ABiM.405 = c("gex_type", "sample", "library_protocol", "ssp_er", "ssp_pr", "ssp_her2", "ssp_her2by_e_rssp", "ssp_ki67", "ssp_nhg"),
    Normal.66 = c("gex_type", "fraction_duplication", "aligned_pairs", "library_protocol", "rna_qc_rin_rqs", "c_clust"),
    OSLO2EMIT0.103 = c("gex_type", "library_protocol", "ssp_cc15", "ssp_pam50", "ssp_subtype", "ssp_ror_as_t0", "ssp_ror_num_tadj", "ssp_ror_risk_cat", "ssp_ror_binary_risk_cat", "ssp_etr", "ssp_etr_dc", "prosigna_subtype", "prosigna_ror", "prosigna_ror_risk_cat", "prosigna_ror_binary_risk_cat", "prosigna_etr", "prosigna_etr_dc"),
    SCANB.9206 = c("gex_type", "follow_up_cohort", "test_set", "training_set", "merged", "library", "rna", "sample", "case", "patient",
    "fraction_duplication", "aligned_pairs", "pm_reads", "pt_reads", "pf_reads", "read_string", "sequencer_serial", "pool_name",
    "library_protocol", "library_barcode", "library_batch_no", "rna_qc_rin_rqs", "rna_nd_conc", "qiacube_batch_no", "reference_year",
    "reference_source", "sampling_date", "days_to_lab", "minutes_to_rna_later", "clin_group", "eval_group", "eval_group_e_rp_her2n_l_nn50", "c_clust", "ssp_cc15", "ncn_pam50", "ncn_subtype", "ncn_ror_as_t0", "ncn_ro_rsd_as_t0", "ncn_ror_as_t1", "ncn_ro_rsd_as_t1", "ncn_ror_risk_cat", "ncn_ror_binary_risk_cat", "ncn_etr", "ncn_etr_dc", "ssp_pam50", "ssp_subtype", "ssp_ror_as_t0", "ssp_ror_risk_cat", "ssp_ror_binary_risk_cat", "ssp_etr", "ssp_etr_dc", "ssp_er", "ssp_pr", "ssp_ki67", "ssp_nhg", "ssp_her2", "ssp_her2_e_rp", "ssp_her2_e_rn", "ssp_her2by_e_rssp")
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

    if (!is.null(data_file$numSummary) && nrow(data_file$numSummary) >= 1) {
        write_tsv(data_file$numSummary, file.path(paste0(meta_summaries_dir, sheet, "_num.tsv")))            
    }

    if (!is.null(data_file$charSummary) && nrow(data_file$charSummary) >= 1) {
        write_tsv(data_file$charSummary, file.path(paste0(meta_summaries_dir, sheet, "_char.tsv")))
    }
}

meta_file_path %>%
    excel_sheets() %>%
    map(clean_metadata, path = meta_file_path)
