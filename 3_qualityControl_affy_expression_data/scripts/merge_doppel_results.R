# script to combine all doppelgang files.
out_file_dir <- "/Data/doppelgang_results"

file_paths <- list.files(out_file_dir, pattern = "*.tsv", full.names = T)

merged_doppel <- read_tsv(file_paths) %>%
  dplyr::select(sample1, sample2, expr.similarity) %>%
  separate(sample1, c("gseID_1", "Sample_ID_1"), sep = ":") %>%
  separate(sample2, c("gseID_2", "Sample_ID_2"), sep = ":") %>%
  mutate(likely_duplicate = ifelse(expr.similarity > 0.99, "Yes", "No")) %>%
  distinct()
  
write_tsv(merged_doppel, "/Data/merged_doppelgang_results.tsv.gz")
