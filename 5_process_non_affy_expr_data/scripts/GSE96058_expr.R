# This is not a microarray study, so no series matrix file with expression values for getGEO to download.
# Thus we have to download the expression data directly

# get expression data
GSE96058 <- getGEOSuppFiles(GEO = "GSE96058", makeDirectory = F, baseDir = tmp_dir, filter_regex = "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")
tmp <- rownames(GSE96058)
GSE96058_expr_table <- read_csv(tmp) %>%
  rename("HGNC_Symbol" = "...1") %>%
  mutate(across(where(is.double), as.character))

# this block of code takes the sample names from the metadada,
# binds them to the expression matrix which does not have any sample names

gseID <- getGEO("GSE96058")
df_HiSeq <- gseID[[1]]
df_NextSeq <- gseID[[2]]

add_sample_names <- function(data_file) {
  metadata <- pData(data_file)
  gsm_id <- metadata[, c(1, 2)] %>%
  t() %>%
  as_tibble %>%
  rownames_to_column(var = "rowname")
  gsm_id[1, 1] <- "HGNC_Symbol"
  gsm_id[2, 1] <- "HGNC_Symbol"

  gsm_id <- gsm_id %>%
  row_to_names(1, remove_row = TRUE)

  select_data <- dplyr::select(GSE96058_expr_table, all_of(colnames(gsm_id)))
  named_data <- bind_rows(gsm_id, select_data) %>%
    dplyr::select(-contains("repl")) %>%
    row_to_names(1, remove_row = TRUE)

  return(named_data)
}

HiSeq <- add_sample_names(df_HiSeq)
NextSeq <- add_sample_names(df_NextSeq)

print("Writing GSE96058 to file!")
write_tsv(HiSeq, paste0(data_dir, "GSE96058_HiSeq.tsv.gz"))
write_tsv(NextSeq, paste0(data_dir, "GSE96058_NextSeq.tsv.gz"))
