

# Download and save file
download.file("https://dcc.icgc.org/api/v1/download?fn=/current/Projects/BRCA-KR/exp_seq.BRCA-KR.tsv.gz",
              destfile = paste0(tmp_dir, "exp_seq_BRCA_KR.tsv.gz"))
exp_array_BRCA_KR <- read_tsv(paste0(tmp_dir, "exp_seq_BRCA_KR.tsv.gz"))


# clean up data
BRCA_KR <- exp_array_BRCA_KR %>%
  rename("HGNC_Symbol" = "gene_id") %>%
  dplyr::filter(HGNC_Symbol != "SLC35E2") %>%
  dplyr::select(HGNC_Symbol, icgc_donor_id, normalized_read_count) %>%
  mutate(normalized_read_count = log2(normalized_read_count + 1)) %>%
  pivot_wider(names_from = icgc_donor_id, values_from = normalized_read_count)

print("Writing BRCA_KR to file!")
write_tsv(BRCA_KR, paste0(data_dir, "ICGC_KR.tsv.gz"))


# This code snippet was used to identify the gene "SLC35E2" as being duplicated
# donor_id <- exp_array_BRCA_KR %>%
#   dplyr::select(icgc_donor_id) %>%
#   distinct() %>%
#   pull(icgc_donor_id)
# This gave us 50 donors

# gene_id <- exp_array_BRCA_KR %>%
#   dplyr::select(gene_id) %>%
#   distinct() %>%
#   pull(gene_id)
# This gave us 26,730 genes

# both <- exp_array_BRCA_KR %>%
#    dplyr::select(gene_id, icgc_donor_id) %>%
#    distinct() %>%
#    nrow()
# we expected to get 1,336,500 (50 * 26,730) but instead we got 1,336,550 . This told us we had an extra gene

# count_test <- exp_array.BRCA_KR %>%
#   dplyr::select(gene_id, icgc_donor_id) %>%
#   group_by(gene_id) %>%
#   summarise(count = n())
# We were able to identify the gene here because it had a count of 100, whereas others had a count of 50
# we then modified the code above and added a filter step
