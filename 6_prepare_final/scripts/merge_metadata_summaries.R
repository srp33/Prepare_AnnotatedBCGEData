
#  The code below summerizes all of the metadata summaries in one file
datadir <- "/Data/metadata_summaries/"
file_paths <- list.files(datadir, full.names = T)
out_file_path <- "/Data/merged_metadata_summary.tsv"

big_df <- NULL

for (file in file_paths) {
    df <- read_tsv(file)

    if (is.null(big_df)) {
        big_df <- df
    } else {
        big_df <- bind_rows(big_df, df)
    }
}
big_df <- big_df[order(big_df$Variable), ]

big_df <- big_df %>%
    distinct(Dataset_ID, Variable, .keep_all = TRUE)

write_tsv(as_tibble(big_df), out_file_path)
print(paste0("Saved to ", out_file_path))


# #  The code below summerizes all of the metadata variables in one file
# datadir <- "/Data/analysis_ready_metadata"
# file_paths <- list.files(datadir, full.names = T)
# out_file_path <- "/Data/merged_metadata.tsv"

# big_column_names <- NULL

# for (file in file_paths) {
#     df <- read_tsv(file)
#     column_names <- as_tibble(colnames(df)) %>%
#         mutate(df[1, 1], .before = value)

#     if (is.null(big_column_names)) {
#         big_column_names <- column_names
#     } else {
#         big_column_names <- bind_rows(big_column_names, column_names)
#     }
# }

# big_column_names <- big_column_names[order(big_column_names$value), ]

# write_tsv(as_tibble(big_column_names), out_file_path)
# print(paste0("Saved to ", out_file_path))
