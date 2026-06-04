# This is not a microarray study, so no series matrix file with expression values for getGEO to download.
# Thus we have to download the expression data directly

# get expression data
GSE81538 <- getGEOSuppFiles(GEO = "GSE81538", makeDirectory = F, baseDir = tmp_dir, filter_regex = "GSE81538_gene_expression_405_transformed.csv.gz")
tmp <- rownames(GSE81538)
GSE81538_expr_table <- read_csv(tmp) %>%
  rename("HGNC_Symbol" = "...1") %>%
  mutate(across(where(is.double), as.character))

# this block of code takes the sample names from the metadata,
# binds them to the expression matrix which doesn't have any sample names

gseID <- getGEO("GSE81538")
df <- gseID[[1]]
metadata <- pData(df)
gsm_id <- metadata[, c(1, 2)] %>%
  t() %>%
  as_tibble %>%
  rownames_to_column(var = "rowname")

gsm_id[1, 1] <- "HGNC_Symbol"
gsm_id[2, 1] <- "HGNC_Symbol"

gsm_id <- gsm_id %>%
  row_to_names(1, remove_row = TRUE)

GSE81538_expr_df <- bind_rows(gsm_id, GSE81538_expr_table) %>%
  row_to_names(1, remove_row = TRUE)

print("Writing GSE81538 to file!")
write_tsv(GSE81538_expr_df, paste0(data_dir, "GSE81538.tsv.gz"))
