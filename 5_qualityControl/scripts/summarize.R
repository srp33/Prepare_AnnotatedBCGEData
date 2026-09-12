library(tidyverse)
library(data.table)

# We don't need to do anything with doppelgangR smoking gun results because it did
#   not identifying any smoking guns.

# fread is much faster than read_tsv for thousands of small files.
# Comment lines (#) from SIS/Jaccard headers are stripped via grep.
doppelgangR_metadata_files <- list.files(
  "/Data/doppelgangR_metadata",
  pattern = "____samples\\.tsv\\.gz$",
  full.names = TRUE
)

# Start with an empty tibble that already has the expected columns so
# bind_rows() does not fail on the first (0-column) bind.
doppelgangR_metadata <- tibble(
  sample_id1 = character(),
  sample_id2 = character(),
  n_shared_values = double(),
  max_possible_shared_values = double(),
  shared_column_values = character()
)

for (file_path in doppelgangR_metadata_files) {
  file_data <- tryCatch(
    fread(
      cmd = sprintf("gzip -cd %s | grep -v '^#'", shQuote(file_path)),
      sep = "\t",
      header = TRUE,
      na.strings = c("", "NA"),
      # Prevent empty shared_column_values from being guessed as logical.
      colClasses = list(
        character = c("sample_id1", "sample_id2", "shared_column_values"),
        numeric = c(
          "shared_information_score",
          "n_shared_values",
          "max_possible_shared_values"
        )
      )
    ) %>% as_tibble(),
    error = function(e) NULL
  )

  # Skip unreadable or header-only files.
  if (is.null(file_data) || nrow(file_data) == 0) {
    next
  }

  file_data <- file_data %>%
    mutate(
      sample_id1 = as.character(sample_id1),
      sample_id2 = as.character(sample_id2),
      shared_column_values = as.character(shared_column_values),
      n_shared_values = as.numeric(n_shared_values),
      max_possible_shared_values = as.numeric(max_possible_shared_values)
    ) %>%
    select(
      sample_id1,
      sample_id2,
      n_shared_values,
      max_possible_shared_values,
      shared_column_values
    )

  doppelgangR_metadata <- bind_rows(doppelgangR_metadata, file_data)
}

# Drop reverse duplicates: (A, B) and (B, A) count as the same pair.
doppelgangR_metadata <- doppelgangR_metadata %>%
  mutate(
    id_a = pmin(sample_id1, sample_id2),
    id_b = pmax(sample_id1, sample_id2)
  ) %>%
  distinct(id_a, id_b, .keep_all = TRUE) %>%
  select(-id_a, -id_b)

print(doppelgangR_metadata, width = Inf)
print(dim(doppelgangR_metadata))

#doppelgangR_expr_data <- read_tsv("/Data/doppelgangR_expr_data/")
