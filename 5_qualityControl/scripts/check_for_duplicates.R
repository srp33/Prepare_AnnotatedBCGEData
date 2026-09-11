# ============================================================
# Shared Information Score (SIS) for cross-dataset sample matching
# ============================================================
#
# The algorithm finds sample pairs that share rare, specific metadata
# values across two datasets, even when column names differ. For each
# pair of samples, it finds values that appear in both samples (in any
# columns). A value that matches in multiple columns counts once per
# column on the smaller side: min(# columns with that value in sample 1,
# # columns with that value in sample 2). So ER=Positive and PR=Positive
# are two matches, not one. Each match is weighted by how rare that
# value is. Matching on a less common value raises the score, because
# two unrelated samples are less likely to share it by chance.
# Example: two samples that both have a rare mutation score higher
# than two that only both say "female", even if each pair matches on
# the same number of fields. A value present in every sample adds
# nothing. These weights are summed into a Shared Information Score.
# Higher scores mean the pair is a more likely duplicate.
# Only pairs with at least 3 shared values, an overlap fraction of at
# least 0.5, and different sample IDs are written to the output file.
# ============================================================

datadir <- "/Data/expression_data4"
metadata_dir <- "/Data/prelim_metadata2"

getExprData <- function(expr_file_path) {
  print(paste0("Reading ", expr_file_path))

  expr_data <- read_tsv(expr_file_path)
  gene_names <- pull(expr_data, Entrez_Gene_ID)
  expr_data <- dplyr::select(expr_data, -Dataset_ID, -Entrez_Gene_ID, -Ensembl_Gene_ID, -HGNC_Symbol, -Chromosomal_Band) %>%
    as.matrix()
  rownames(expr_data) <- gene_names

  return(expr_data)
}

getMetadata <- function(metadata_file_path) {
  print(paste0("Reading ", metadata_file_path))

  metadata <- read_tsv(metadata_file_path) %>%
    dplyr::select(-Dataset_ID, -Platform_ID)

  if (ncol(metadata) == 0) {
    return(NULL)
  } else {
    sampleIDs <- dplyr::pull(metadata, Sample_ID)
    metadata <- dplyr::select(metadata, -Sample_ID) %>%
      as.data.frame()
    rownames(metadata) <- sampleIDs

    return(metadata)
  }
}

convertIntStringToFloatString <- function(x) {
  is_integer_string <- grepl("^-?[0-9]+$", x)

  x[is_integer_string] <- sprintf("%.1f", as.numeric(x[is_integer_string]))

  return(x)
}

# Normalize a metadata data frame to a character matrix so that
# "45" and "45.0" are treated as the same value.
normalize_metadata_values <- function(metadata) {
  mat <- as.matrix(metadata)
  storage.mode(mat) <- "character"
  mat[is.na(metadata)] <- NA_character_
  mat[] <- convertIntStringToFloatString(mat)
  mat
}

# Build a named list: sample_id -> unique non-missing metadata values.
sample_value_sets <- function(value_mat) {
  setNames(
    lapply(seq_len(nrow(value_mat)), function(i) {
      unique(na.omit(value_mat[i, ]))
    }),
    rownames(value_mat)
  )
}

# For each sample, map each value to the column name(s) where it appears.
sample_value_to_columns <- function(value_mat) {
  cols <- colnames(value_mat)
  setNames(
    lapply(seq_len(nrow(value_mat)), function(i) {
      row <- value_mat[i, ]
      ok <- !is.na(row)
      if (!any(ok)) {
        return(list())
      }
      split(unname(cols[ok]), row[ok])
    }),
    rownames(value_mat)
  )
}

# Count how many samples (across both datasets) contain each value in
# at least one column. Rarity weight is -log2(proportion of samples),
# so ubiquitous values get weight ~0 and rare values get large weight.
value_rarity_weights <- function(value_sets1, value_sets2) {
  n_total <- length(value_sets1) + length(value_sets2)
  value_counts <- table(unlist(c(value_sets1, value_sets2), use.names = FALSE))
  # IDF-style weight; values present in every sample contribute 0.
  setNames(-log2(as.numeric(value_counts) / n_total), names(value_counts))
}

# Format one shared value for the detail column:
#   col1a,col1b=value|col2a,col2b
format_shared_match <- function(cols1, cols2, value) {
  paste0(
    paste(cols1, collapse = ","),
    "=", value, "|",
    paste(cols2, collapse = ",")
  )
}

# ------------------------------------------------------------
# Shared Information Score for all sample pairs between two datasets.
#
# Uses an inverted index over values: for each unique value, add its
# rarity weight to every sample pair that both carry that value. This
# is equivalent to set-intersection scoring but avoids an O(n1*n2)
# pass over empty pairs when most values are sparse.
# ------------------------------------------------------------
calcSharedInformationScores <- function(metadata1, metadata2) {
  value_mat1 <- normalize_metadata_values(metadata1)
  value_mat2 <- normalize_metadata_values(metadata2)

  value_sets1 <- sample_value_sets(value_mat1)
  value_sets2 <- sample_value_sets(value_mat2)
  value_cols1 <- sample_value_to_columns(value_mat1)
  value_cols2 <- sample_value_to_columns(value_mat2)
  weights <- value_rarity_weights(value_sets1, value_sets2)

  # Invert: value -> sample IDs that contain it in each dataset.
  samples_by_value1 <- list()
  for (sid in names(value_sets1)) {
    for (v in value_sets1[[sid]]) {
      samples_by_value1[[v]] <- c(samples_by_value1[[v]], sid)
    }
  }
  samples_by_value2 <- list()
  for (sid in names(value_sets2)) {
    for (v in value_sets2[[sid]]) {
      samples_by_value2[[v]] <- c(samples_by_value2[[v]], sid)
    }
  }

  # Only values present in both datasets can contribute to any pair.
  shared_values <- intersect(names(samples_by_value1), names(samples_by_value2))

  score_matrix <- matrix(0, nrow = nrow(metadata1), ncol = nrow(metadata2))
  rownames(score_matrix) <- rownames(metadata1)
  colnames(score_matrix) <- rownames(metadata2)

  match_count_matrix <- matrix(0, nrow = nrow(metadata1), ncol = nrow(metadata2))
  rownames(match_count_matrix) <- rownames(metadata1)
  colnames(match_count_matrix) <- rownames(metadata2)

  # Character matrix of shared column/value match descriptions.
  detail_matrix <- matrix("", nrow = nrow(metadata1), ncol = nrow(metadata2))
  rownames(detail_matrix) <- rownames(metadata1)
  colnames(detail_matrix) <- rownames(metadata2)

  for (v in shared_values) {
    w <- unname(weights[v])
    if (is.na(w)) next

    s1 <- samples_by_value1[[v]]
    s2 <- samples_by_value2[[v]]
    # Count column-level matches, not unique strings. ER=Positive and
    # PR=Positive are two matches when both sides have that value in
    # two columns: min(n_cols1, n_cols2).
    n_cols1 <- vapply(s1, function(a) length(value_cols1[[a]][[v]]), integer(1))
    n_cols2 <- vapply(s2, function(b) length(value_cols2[[b]][[v]]), integer(1))
    n_match <- outer(n_cols1, n_cols2, pmin)
    # Ubiquitous values have weight 0 and do not change the score,
    # but still count toward n_shared_values and appear in the detail column.
    if (w > 0) {
      score_matrix[s1, s2] <- score_matrix[s1, s2] + w * n_match
    }
    match_count_matrix[s1, s2] <- match_count_matrix[s1, s2] + n_match

    for (a in s1) {
      cols_a <- value_cols1[[a]][[v]]
      for (b in s2) {
        piece <- format_shared_match(cols_a, value_cols2[[b]][[v]], v)
        if (detail_matrix[a, b] == "") {
          detail_matrix[a, b] <- piece
        } else {
          detail_matrix[a, b] <- paste(detail_matrix[a, b], piece, sep = ";")
        }
      }
    }
  }

  score_df <- as.data.frame(as.table(score_matrix))
  colnames(score_df) <- c("sample_id1", "sample_id2", "shared_information_score")
  match_df <- as.data.frame(as.table(match_count_matrix))
  colnames(match_df) <- c("sample_id1", "sample_id2", "n_shared_values")
  # as.table() is numeric-only; expand the character detail matrix manually.
  detail_df <- data.frame(
    sample_id1 = rep(rownames(detail_matrix), times = ncol(detail_matrix)),
    sample_id2 = rep(colnames(detail_matrix), each = nrow(detail_matrix)),
    shared_column_values = as.vector(detail_matrix),
    stringsAsFactors = FALSE
  )

  # Ceiling is the smaller number of non-missing fields, since each
  # match consumes one column on each side.
  n_fields1 <- setNames(rowSums(!is.na(value_mat1)), rownames(value_mat1))
  n_fields2 <- setNames(rowSums(!is.na(value_mat2)), rownames(value_mat2))

  score_df %>%
    inner_join(match_df, by = c("sample_id1", "sample_id2")) %>%
    inner_join(detail_df, by = c("sample_id1", "sample_id2")) %>%
    mutate(
      max_possible_shared_values = pmin(
        n_fields1[as.character(sample_id1)],
        n_fields2[as.character(sample_id2)]
      ),
      shared_information_score = round(shared_information_score, 6)
    ) %>%
    relocate(max_possible_shared_values, .after = n_shared_values) %>%
    # Keep only pairs with enough shared information to be interesting.
    filter(
      n_shared_values >= 3,
      max_possible_shared_values > 0,
      (n_shared_values / max_possible_shared_values) >= 0.5,
      as.character(sample_id1) != as.character(sample_id2)
    ) %>%
    arrange(desc(shared_information_score), sample_id1, sample_id2)
}

sis_output_comment <- c(
  "# Shared Information Score (SIS) for candidate duplicate samples.",
  "# Higher scores = stronger evidence the two samples are the same individual.",
  "# Score = sum of (rarity weight * column-level matches) for shared values.",
  "# A less common shared value raises the score, because unrelated samples",
  "# are less likely to match on it by chance. Example: sharing a rare mutation",
  "# scores higher than both saying 'female', even with the same number of matches.",
  "# A value in k columns on each side counts k times (e.g. ER+ and PR+ are two matches).",
  "# n_shared_values is the number of column-level matches (sum of min(#cols1, #cols2) per value).",
  "# max_possible_shared_values = min(# non-missing fields in sample1, # in sample2).",
  "# shared_column_values lists matches as: col1[,col1b]=value|col2[,col2b] (semicolon-separated).",
  "# Only pairs with n_shared_values >= 3, overlap fraction >= 0.5, and sample_id1 != sample_id2 are saved."
)

jaccard_output_comment <- c(
  "# Metadata column pairs across two datasets (col1 from dataset 1, col2 from dataset 2).",
  "# jaccard_score measures how similar the value *distributions* are (not just shared labels):",
  "#   sum(min(p1,p2)) / sum(max(p1,p2)); 1 = same frequencies, 0 = no shared values.",
  "# High scores suggest the two columns may represent the same kind of variable.",
  "# This file is for inspection only; SIS sample matching does not filter on these scores."
)

# write_tsv/vroom requires a binary connection, so for commented .tsv.gz
# output we write via a text gzfile connection + write.table instead.
write_commented_tsv_gz <- function(df, path, comments) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  writeLines(comments, con)
  write.table(
    df,
    file = con,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = ""
  )
}

write_sis_samples <- function(df, path) {
  write_commented_tsv_gz(df, path, sis_output_comment)
}

write_jaccard_variables <- function(df, path) {
  write_commented_tsv_gz(df, path, jaccard_output_comment)
}

empty_sis_tbl <- function() {
  tibble(
    sample_id1 = character(),
    sample_id2 = character(),
    shared_information_score = numeric(),
    n_shared_values = integer(),
    max_possible_shared_values = integer(),
    shared_column_values = character()
  )
}

# Column-pair similarity based on value *distributions*, not just the
# set of unique labels. Unique-set Jaccard scores 1 whenever two columns
# share the same labels (e.g. both {1,2,3} or both {yes,no}), even if one
# is nearly all "1" and the other is uniform — which incorrectly marks
# unrelated columns as matches. Here we compare relative frequencies
# (Ruzicka / probability Jaccard): sum(min(p1,p2)) / sum(max(p1,p2)).
calcJaccardScore <- function(metadata1, metadata2, col1_vector, col2_vector) {
  vapply(seq_along(col1_vector), function(i) {
    x <- convertIntStringToFloatString(as.character(metadata1[[col1_vector[i]]]))
    y <- convertIntStringToFloatString(as.character(metadata2[[col2_vector[i]]]))
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]

    if (length(x) == 0 || length(y) == 0) {
      return(NA_real_)
    }

    px <- table(x)
    py <- table(y)
    px <- px / sum(px)
    py <- py / sum(py)

    all_vals <- union(names(px), names(py))
    p1 <- setNames(rep(0, length(all_vals)), all_vals)
    p2 <- p1
    p1[names(px)] <- as.numeric(px)
    p2[names(py)] <- as.numeric(py)

    sum(pmin(p1, p2)) / sum(pmax(p1, p2))
  }, numeric(1))
}

processCombo <- function(file_path1, file_path2, dataset_id1, dataset_id2, metadata_file_path1, metadata_file_path2, sg_out_file_path, md_out_file_path, ed_out_file_path) {
  if (file.exists(md_out_file_path)) {
    return(NULL)
  }

  expr_data1 <- getExprData(file_path1)
  expr_data2 <- getExprData(file_path2)

  metadata1 <- getMetadata(metadata_file_path1)
  metadata2 <- getMetadata(metadata_file_path2)

  sample_ids1 <- sort(intersect(colnames(expr_data1), rownames(metadata1)))
  sample_ids2 <- sort(intersect(colnames(expr_data2), rownames(metadata2)))

  if (length(sample_ids1) < 5) {
    stop(paste0("There are few, if any, matching samples between the metadata and expression data for ", file_path1, "."))
  }
  if (length(sample_ids2) < 5) {
    stop(paste0("There are few, if any, matching samples between the metadata and expression data for ", file_path2, "."))
  }

  genes <- sort(intersect(rownames(expr_data1), rownames(expr_data2)))

  if (length(genes) < 1000) {
    stop(paste0("The number of genes overlapping between ", dataset_id1, " and ", dataset_id2, " is less than 1000 [", length(genes), "], so there must be a problem."))
  }

  expr_data1 <- expr_data1[genes,sample_ids1]
  expr_data2 <- expr_data2[genes,sample_ids2]

  metadata1 <- metadata1[sample_ids1, , drop = FALSE]
  metadata2 <- metadata2[sample_ids2, , drop = FALSE]

  if (!file.exists(sg_out_file_path)) {
    # We use only a few genes to reduce memory usage.
    # This is fine because we don't use the expression data here,
    # but doppelgangR requires it to be in the ExpressionSet.
    eSet1 <- ExpressionSet(assayData = expr_data1[1:5,], phenoData = AnnotatedDataFrame(metadata1))
    eSet2 <- ExpressionSet(assayData = expr_data2[1:5,], phenoData = AnnotatedDataFrame(metadata2))

    dopple_data <- list(eSet1, eSet2)
    names(dopple_data) <- c(dataset_id1, dataset_id2)

    result <- doppelgangR(
      dopple_data,
      corFinder.args = NULL,
      phenoFinder.args = NULL,
      automatic.smokingguns = TRUE
    )

    result <- summary(result)
    result <- dplyr::select(result, sample1, sample2, smokinggun.similarity, smokinggun.doppel)
    write_tsv(result, sg_out_file_path)
  }

  if (!file.exists(md_out_file_path)) {
    if (is.null(metadata1) || is.null(metadata2) ||
        ncol(metadata1) == 0 || ncol(metadata2) == 0) {
      write_sis_samples(empty_sis_tbl(), md_out_file_path)
    } else {
      print(paste0("Calculating Shared Information Scores for ", dataset_id1, " and ", dataset_id2))
      write_sis_samples(
        calcSharedInformationScores(metadata1, metadata2),
        md_out_file_path
      )

      # Variable-pair Jaccard scores are written separately for inspection;
      # SIS itself compares values across all column combinations.
      candidate_metadata_combos <- expand.grid(
        col1 = colnames(metadata1),
        col2 = colnames(metadata2),
        stringsAsFactors = FALSE
      )
      if (nrow(candidate_metadata_combos) > 0) {
        candidate_metadata_combos <- mutate(
          candidate_metadata_combos,
          jaccard_score = calcJaccardScore(metadata1, metadata2, col1, col2)
        ) %>%
          arrange(desc(jaccard_score), col1, col2)
        write_jaccard_variables(
          candidate_metadata_combos,
          sub("____samples.tsv.gz", "____variables.tsv.gz", md_out_file_path)
        )
      }
    }
  }

  if (!file.exists(ed_out_file_path)) {
    print(paste0("Calculating expression correlation for ", dataset_id1, " and ", dataset_id2))
    expr_data <- cbind(expr_data1, expr_data2)

    cor_matrix <- cor(expr_data, method = "spearman")

    cor_tbl <- cor_matrix |>
      as.data.frame() |>
      rownames_to_column("sample1") |>
      mutate(row_num = row_number()) |>
      pivot_longer(
        -c(sample1, row_num),
        names_to = "sample2",
        values_to = "correlation_coefficient"
      ) |>
      mutate(col_num = match(sample2, colnames(cor_matrix))) |>
      filter(row_num > col_num) |>
      select(sample1, sample2, correlation_coefficient) |>
      arrange(desc(correlation_coefficient), sample1, sample2)

    write_tsv(cor_tbl, ed_out_file_path)
  }
}

# Enable parallelization
num_parallel = 16
registerDoParallel(num_parallel)
stopifnot(foreach::getDoParWorkers() == num_parallel)

# Run dopplegangR for pairwise comparisons of datasets.
file_paths <- list.files(datadir, full.names = T)

# Get all unique pairs.
pairs <- as.data.frame(
  t(combn(file_paths, 2)),
  stringsAsFactors = FALSE
)

set.seed(0)

colnames(pairs) <- c("file_path1", "file_path2")
pairs <- mutate(pairs, dataset_id1 = basename(file_path1)) %>%
  mutate(dataset_id2 = basename(file_path2)) %>%
  mutate(dataset_id1 = str_replace(dataset_id1, "\\.tsv\\.gz", "")) %>%
  mutate(dataset_id2 = str_replace(dataset_id2, "\\.tsv\\.gz", "")) %>%
  mutate(metadata_file_path1 = str_c(metadata_dir, "/", dataset_id1, ".tsv")) %>%
  mutate(metadata_file_path2 = str_c(metadata_dir, "/", dataset_id2, ".tsv")) %>%
  mutate(sg_out_file_path = str_c("/Data/doppelgangR_smokinggun/", dataset_id1, "_", dataset_id2, ".tsv.gz")) %>%
  mutate(md_out_file_path = str_c("/Data/doppelgangR_metadata/", dataset_id1, "_", dataset_id2, "____samples.tsv.gz")) %>%
  mutate(ed_out_file_path = str_c("/Data/doppelgangR_expr_data/", dataset_id1, "_", dataset_id2, ".tsv.gz")) %>%
  filter(dataset_id1 == "ABiM.100")
#filter(dataset_id1 == "GSE12276" & dataset_id2 == "GSE12763")
# filter(dataset_id1 == "ABiM.100" & dataset_id2 == "ABiM.405")
#filter(dataset_id1 == "SCANB.9206" | dataset_id2 == "SCANB.9206")
#filter(dataset_id1 == "GSE96058_HiSeq" & dataset_id2 == "SCANB.9206")
 # slice_sample(prop = 1)

foreach (i = 1:nrow(pairs)) %dopar% {
#for (i in 1:nrow(pairs)) {
  row <- as.vector(as.matrix(pairs[i,]))
  file_path1 <- row[1]
  file_path2 <- row[2]
  dataset_id1 <- row[3]
  dataset_id2 <- row[4]
  metadata_file_path1 <- row[5]
  metadata_file_path2 <- row[6]
  sg_out_file_path <- row[7]
  md_out_file_path <- row[8]
  ed_out_file_path <- row[9]

  processCombo(file_path1, file_path2, dataset_id1, dataset_id2, metadata_file_path1, metadata_file_path2, sg_out_file_path, md_out_file_path, ed_out_file_path)
}

unlink("cache", recursive = TRUE, force = TRUE)
