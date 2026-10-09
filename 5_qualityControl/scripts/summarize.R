library(tidyverse)
library(data.table)

options(readr.show_col_types = FALSE, readr.show_progress = FALSE)

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

# Gene-level values for one sample pair, so identical and differing genes can be inspected.
write_pair_expression <- function(dataset_id1, dataset_id2, sample1, sample2) {
  data1 <- read_tsv(
    paste0("/Data/expression_data4/", dataset_id1, ".tsv.gz"),
    col_select = all_of(c("Entrez_Gene_ID", sample1))
  )
  data2 <- read_tsv(
    paste0("/Data/expression_data4/", dataset_id2, ".tsv.gz"),
    col_select = all_of(c("Entrez_Gene_ID", sample2))
  )

  genes <- intersect(data1$Entrez_Gene_ID, data2$Entrez_Gene_ID)
  data1 <- data1[match(genes, data1$Entrez_Gene_ID), , drop = FALSE]
  data2 <- data2[match(genes, data2$Entrez_Gene_ID), , drop = FALSE]

  value1 <- data1[[sample1]]
  value2 <- data2[[sample2]]
  pair_values <- tibble(
    Entrez_Gene_ID = genes,
    !!sample1 := value1,
    !!sample2 := value2,
    same = round(value1, 3) == round(value2, 3)
  ) %>%
    # Unequal values, including missing comparisons, come before matches.
    arrange(coalesce(same, FALSE), Entrez_Gene_ID)

  out_path <- paste0(
    "/Data/", dataset_id1, "_", sample1, "__", dataset_id2, "_", sample2, ".tsv"
  )
  write_tsv(pair_values, out_path)

  plot_path <- sub("\\.tsv$", ".pdf", out_path)
  scatter <- ggplot(pair_values, aes(x = .data[[sample1]], y = .data[[sample2]])) +
    geom_point(alpha = 0.4, size = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(
      x = paste(dataset_id1, sample1),
      y = paste(dataset_id2, sample2)
    ) +
    theme_bw()
  ggsave(plot_path, scatter, width = 6, height = 6)

  c(tsv = out_path, pdf = plot_path)
}

# doppelgangR_metadata <- get_doppelgangR_metadata()
expr_pairs <- get_doppelgangR_expr_data()

# Read each expression file once, and only the sample columns used in pairs.
needed_samples <- bind_rows(
  expr_pairs %>% transmute(dataset_id = dataset_id1, sample_id = sample1),
  expr_pairs %>% transmute(dataset_id = dataset_id2, sample_id = sample2)
) %>%
  distinct() %>%
  summarise(sample_ids = list(sample_id), .by = dataset_id)

expr_by_dataset <- set_names(
  lapply(seq_len(nrow(needed_samples)), function(i) {
    dataset_id <- needed_samples$dataset_id[i]
    print(paste0("Reading expression data for ", dataset_id))
    read_tsv(
      paste0("/Data/expression_data4/", dataset_id, ".tsv.gz"),
      col_select = all_of(c("Entrez_Gene_ID", needed_samples$sample_ids[[i]]))
    )
  }),
  needed_samples$dataset_id
)

dataset_pairs <- expr_pairs %>%
  distinct(dataset_id1, dataset_id2)

results <- vector("list", nrow(dataset_pairs))
for (i in seq_len(nrow(dataset_pairs))) {
  dataset_id1 <- dataset_pairs$dataset_id1[i]
  dataset_id2 <- dataset_pairs$dataset_id2[i]
  print(paste0(
    "Comparing ", dataset_id1, " and ", dataset_id2,
    " (", i, " of ", nrow(dataset_pairs), ")"
  ))

  sample_pairs <- expr_pairs %>%
    filter(dataset_id1 == .env$dataset_id1, dataset_id2 == .env$dataset_id2)

  # Align genes once, then compare every sample pair from these two datasets.
  data1 <- expr_by_dataset[[dataset_id1]]
  data2 <- expr_by_dataset[[dataset_id2]]
  genes <- intersect(data1$Entrez_Gene_ID, data2$Entrez_Gene_ID)
  data1 <- data1[match(genes, data1$Entrez_Gene_ID), , drop = FALSE]
  data2 <- data2[match(genes, data2$Entrez_Gene_ID), , drop = FALSE]

  # NA is neither same nor different. Equal values count as the same.
  counts <- vapply(seq_len(nrow(sample_pairs)), function(j) {
    x <- round(data1[[sample_pairs$sample1[j]]], 3)
    y <- round(data2[[sample_pairs$sample2[j]]], 3)
    c(
      num_same = as.integer(sum(x == y, na.rm = TRUE)),
      num_different = as.integer(sum(x != y, na.rm = TRUE))
    )
  }, integer(2))

  results[[i]] <- sample_pairs %>%
    mutate(
      num_same = counts["num_same", ],
      num_different = counts["num_different", ]
    )
}

results <- bind_rows(results) %>%
  arrange(desc(correlation_coefficient), num_different)

write_tsv(results, "/Data/test.tsv")

# write_pair_expression("ABiM.100", "ABiM.405", "ABiM100.001", "ABiM405.269")
# write_pair_expression("E_TABM_158", "GSE7378", "b0341", "GSM177501")
# write_pair_expression("GSE20194", "GSE25055", "GSM505578", "GSM615351")
# write_pair_expression("GSE20194", "GSE25055", "GSM505565", "GSM615337")
# write_pair_expression("GSE21653", "GSE31448", "GSM540175", "GSM781316")
write_pair_expression("GSE23720", "GSE31448", "GSM585333", "GSM781316")
write_pair_expression("GSE96058_HiSeq", "GSE96058_NextSeq", "GSM2530337", "GSM2530057")


#TODO: Filter sample pairs based on this. Identify criteria for filtering based on one or the other, not necessarily both.

#TODO: Create code similar to what is above for aggregating the Variables.tsv.gz files. Then filter metadata variables based on the values in the Variables.tsv.gz files.
#      Need to have a way of picking one when two variables are perfectly overlapping. (First make sure the logic right.)

#TODO: Save files that indicate what needs to be filtered out. Use that in a future step to filter the data and ontology mappings.

#TODO: Remove the IQRay stuff from the Dockerfile and from the pipeline.
