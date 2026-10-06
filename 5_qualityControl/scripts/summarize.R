library(tidyverse)
library(data.table)

# We don't need to do anything with doppelgangR smoking gun results because it did
#   not identifying any smoking guns.
#
# Pair files written by check_for_duplicates.R include dataset_id1 / dataset_id2
# columns (sample1/sample_id1 from dataset_id1, sample2/sample_id2 from dataset_id2).

# Stop immediately on a truncated or otherwise invalid .gz file.
# fread can otherwise keep going after gzip writes "unexpected end of file".
assert_gzip_ok <- function(file_path) {
  status <- system2("gzip", c("-t", file_path))
  if (!identical(status, 0L)) {
    stop("Cannot read ", file_path, call. = FALSE)
  }
}

get_doppelgangR_metadata <- function() {
  summary_file_path <- "/Data/doppelgangR_metadata/summary.tsv.gz"

  if (file.exists(summary_file_path)) {
    print(paste0("Reading summary file from ", summary_file_path))
    return(read_tsv(summary_file_path))
  }

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
    dataset_id1 = character(),
    dataset_id2 = character(),
    sample_id1 = character(),
    sample_id2 = character(),
    n_shared_values = double(),
    max_possible_shared_values = double(),
    shared_column_values = character()
  )

  for (file_path in doppelgangR_metadata_files) {
    print(file_path)
    assert_gzip_ok(file_path)
    file_data <- fread(
      cmd = sprintf("gzip -cd %s | grep -v '^#'", shQuote(file_path)),
      sep = "\t",
      header = TRUE,
      na.strings = c("", "NA"),
      # Prevent empty shared_column_values from being guessed as logical.
      colClasses = list(
        character = c(
          "dataset_id1",
          "dataset_id2",
          "sample_id1",
          "sample_id2",
          "shared_column_values"
        ),
        numeric = c(
          "shared_information_score",
          "n_shared_values",
          "max_possible_shared_values"
        )
      )
    ) %>% as_tibble()

    # A readable file can still have no data rows.
    if (nrow(file_data) == 0) {
      next
    }
    if (!all(c("dataset_id1", "dataset_id2") %in% names(file_data))) {
      stop("Missing dataset_id columns in ", file_path, call. = FALSE)
    }

    file_data <- file_data %>%
      mutate(
        dataset_id1 = as.character(dataset_id1),
        dataset_id2 = as.character(dataset_id2),
        sample_id1 = as.character(sample_id1),
        sample_id2 = as.character(sample_id2),
        shared_column_values = as.character(shared_column_values),
        n_shared_values = as.numeric(n_shared_values),
        max_possible_shared_values = as.numeric(max_possible_shared_values)
      ) %>%
      dplyr::select(
        dataset_id1,
        dataset_id2,
        sample_id1,
        sample_id2,
        n_shared_values,
        max_possible_shared_values,
        shared_column_values
      )

    doppelgangR_metadata <- bind_rows(doppelgangR_metadata, file_data)
  }

  # Drop reverse duplicates. A pair is the same if the two
  # (dataset, sample) sides are swapped.
  doppelgangR_metadata <- doppelgangR_metadata %>%
    mutate(
      side_a = paste(dataset_id1, sample_id1, sep = "\t"),
      side_b = paste(dataset_id2, sample_id2, sep = "\t"),
      id_a = pmin(side_a, side_b),
      id_b = pmax(side_a, side_b)
    ) %>%
    distinct(id_a, id_b, .keep_all = TRUE) %>%
    dplyr::select(-side_a, -side_b, -id_a, -id_b)

  write_tsv(doppelgangR_metadata, summary_file_path)

  return(doppelgangR_metadata)
}

get_doppelgangR_expr_data <- function() {
  summary_file_path <- "/Data/doppelgangR_expr_data/summary.tsv.gz"

  if (file.exists(summary_file_path)) {
    print(paste0("Reading summary file from ", summary_file_path))
    return(read_tsv(summary_file_path))
  }

  # fread is much faster than read_tsv for thousands of small files.
  doppelgangR_expr_data_files <- list.files(
    "/Data/doppelgangR_expr_data",
    pattern = "\\.tsv\\.gz$",
    full.names = TRUE
  )
  doppelgangR_expr_data_files <- doppelgangR_expr_data_files[
    basename(doppelgangR_expr_data_files) != "summary.tsv.gz"
  ]

  # Start with an empty tibble that already has the expected columns so
  # bind_rows() does not fail on the first (0-column) bind.
  doppelgangR_expr_data <- tibble(
    dataset_id1 = character(),
    dataset_id2 = character(),
    sample1 = character(),
    sample2 = character(),
    correlation_coefficient = double()
  )

  for (file_path in doppelgangR_expr_data_files) {
    print(file_path)
    assert_gzip_ok(file_path)
    file_data <- fread(
      file_path,
      sep = "\t",
      header = TRUE,
      na.strings = c("", "NA"),
      colClasses = list(
        character = c("dataset_id1", "dataset_id2", "sample1", "sample2"),
        numeric = "correlation_coefficient"
      )
    ) %>% as_tibble()

    # A readable file can still have no data rows.
    if (nrow(file_data) == 0) {
      next
    }
    if (!all(c("dataset_id1", "dataset_id2") %in% names(file_data))) {
      stop("Missing dataset_id columns in ", file_path, call. = FALSE)
    }

    file_data <- file_data %>%
      mutate(
        dataset_id1 = as.character(dataset_id1),
        dataset_id2 = as.character(dataset_id2),
        sample1 = as.character(sample1),
        sample2 = as.character(sample2),
        correlation_coefficient = as.numeric(correlation_coefficient)
      ) %>%
      dplyr::select(dataset_id1, dataset_id2, sample1, sample2, correlation_coefficient)

    doppelgangR_expr_data <- bind_rows(doppelgangR_expr_data, file_data)
  }

  # Drop reverse duplicates. A pair is the same if the two
  # (dataset, sample) sides are swapped.
  doppelgangR_expr_data <- doppelgangR_expr_data %>%
    mutate(
      side_a = paste(dataset_id1, sample1, sep = "\t"),
      side_b = paste(dataset_id2, sample2, sep = "\t"),
      id_a = pmin(side_a, side_b),
      id_b = pmax(side_a, side_b)
    ) %>%
    distinct(id_a, id_b, .keep_all = TRUE) %>%
    dplyr::select(-side_a, -side_b, -id_a, -id_b) %>%
    arrange(desc(correlation_coefficient))

  write_tsv(doppelgangR_expr_data, summary_file_path)

  return(doppelgangR_expr_data)
}

doppelgangR_metadata <- get_doppelgangR_metadata()
# print(doppelgangR_metadata)
print(dim(doppelgangR_metadata))

doppelgangR_expr_data <- get_doppelgangR_expr_data()
# print(doppelgangR_expr_data, n = 100, width = Inf)
print(dim(doppelgangR_expr_data))

#TODO: Filter sample pairs based on this. Identify criteria for filtering based on one or the other, not necessarily both.

#TODO: Create code similar to what is above for aggregating the Variables.tsv.gz files. Then filter metadata variables based on the values in the Variables.tsv.gz files.
#      Need to have a way of picking one when two variables are perfectly overlapping. (First make sure the logic right.)
#TODO: Save files that indicate what needs to be filtered out. Use that in a future step to filter the data and ontology mappings.
#TODO: Remove the IQRay stuff from the Dockerfile and from the pipeline.
