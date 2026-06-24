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

calcJaccardScore <- function(metadata1, metadata2, col1_vector, col2_vector) {
  scores <- c()

  for (i in 1:length(col1_vector)) {
    col1 <- col1_vector[i]
    col2 <- col2_vector[i]

    x <- pull(metadata1, col1)
    y <- pull(metadata2, col2)

    x <- sort(unique(na.omit(as.character(x))))
    y <- sort(unique(na.omit(as.character(y))))

    x <- convertIntStringToFloatString(x)
    y <- convertIntStringToFloatString(y)

    scores <- c(scores, length(intersect(x, y)) / length(union(x, y)))
  }

  return(scores)
}

calcSamplePairScores <- function(dataset_id1, dataset_id2, metadata1, metadata2, metadata_combos) {
  count_matrix <- matrix(0, nrow = nrow(metadata1), ncol = nrow(metadata2))
  rownames(count_matrix) <- rownames(metadata1)
  colnames(count_matrix) <- rownames(metadata2)

  for (i in 1:nrow(metadata_combos)) {
    print(c(dataset_id1, dataset_id2))
    print(metadata_combos[i,])
    col1 <- as.vector(metadata_combos[i,1])
    col2 <- as.vector(metadata_combos[i,2])

    for (j in 1:nrow(metadata1)) {
      sample_id1 <- rownames(metadata1)[j]

      for (k in 1:nrow(metadata2)) {
        sample_id2 <- rownames(metadata2)[k]

        value1 <- metadata1[j,col1]
        value2 <- metadata2[k,col2]

        if (is.null(value1) || is.null(value2) || is.na(value1) || is.na(value2)) {
          next
        }

        if (value1 == value2) {
          count_matrix[sample_id1, sample_id2] <- count_matrix[sample_id1, sample_id2] + 1
        }
      }
    }
  }

  df <- as.data.frame(as.table(count_matrix))
  colnames(df) <- c("sample_id1", "sample_id2", "count")

  filter(df, count > 0) %>%
    arrange(desc(count), sample_id1, sample_id2) %>%
    return()
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
    candidate_metadata_combos <- expand.grid(
      col1 = colnames(metadata1),
      col2 = colnames(metadata2)
    )

    if (nrow(candidate_metadata_combos) == 0) {
      write_tsv(data.frame(sample_id1 = character(), sample_id2 = character(), count = character()), md_out_file_path)
    } else {
      candidate_metadata_combos = mutate(candidate_metadata_combos, jaccard_score = calcJaccardScore(metadata1, metadata2, col1, col2)) %>%
        filter(jaccard_score > 0.1) %>% # This threshold is arbitrary but fairly low by design.
        dplyr::select(-jaccard_score)

      if (nrow(candidate_metadata_combos) == 0) {
        write_tsv(data.frame(sample_id1 = character(), sample_id2 = character(), count = character()), md_out_file_path)
      } else {
        write_tsv(calcSamplePairScores(dataset_id1, dataset_id2, metadata1, metadata2, candidate_metadata_combos), md_out_file_path)
        write_tsv(candidate_metadata_combos, sub("____samples.tsv.gz", "____variables.tsv.gz", md_out_file_path))
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
  mutate(ed_out_file_path = str_c("/Data/doppelgangR_expr_data/", dataset_id1, "_", dataset_id2, ".tsv.gz")) #%>%
#filter(dataset_id1 == "GSE12276" & dataset_id2 == "GSE12763")
#filter(dataset_id1 == "ABiM.100" & dataset_id2 == "ABiM.405")
#filter(dataset_id1 == "SCANB.9206" | dataset_id2 == "SCANB.9206")
#filter(dataset_id1 == "GSE96058_HiSeq" & dataset_id2 == "SCANB.9206")
 # slice_sample(prop = 1)

#foreach (i = 1:nrow(pairs)) %dopar% {
for (i in 1:nrow(pairs)) {
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
