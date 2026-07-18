datadir <- "/Data/expression_data4"
metadata_dir <- "/Data/prelim_metadata2"
out_dir <- "/Data/variation_results"

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
    dplyr::select(-Dataset_ID)

  platform_id <- pull(metadata, Platform_ID) %>%
    unique()

  metadata <- dplyr::select(metadata, -Platform_ID)

  if (ncol(metadata) <= 1) {
    return(list(Platform_ID = platform_id, Metadata = NULL))
  } else {
    sampleIDs <- dplyr::pull(metadata, Sample_ID)
    metadata <- dplyr::select(metadata, -Sample_ID) %>%
      as.data.frame()
    rownames(metadata) <- sampleIDs

    return(list(Platform_ID = platform_id, Metadata = metadata))
  }
}

# ============================================================
# Variance Partitioning QC for Gene Expression Data
# ============================================================
#
# Quantifies how much expression variance each metadata variable
# explains, using a single joint model so correlated variables
# share credit rather than double-counting. The "Unexplained
# variance" row is the fraction not accounted for by any variable
# in the model - a combination of biological noise and any
# unmeasured technical variation such as unrecorded batch effects.
#
# Variance partitioning is fit on a random subset of genes (default
# 1000) after excluding genes with zero variance across samples.
# Canonical correlations between metadata variables are computed first.
# Pairs with CCA = 1 are treated as redundant: the alphabetically second
# variable in each pair is excluded from the variance partitioning model.
#
# Variable roles are assigned automatically:
#   Continuous variables           -> fixed effect (linear)
#   Categorical variables          -> (1|var) random effect
#   Constant columns               -> excluded
#   Near-unique per sample         -> excluded (e.g. SampleID)
#
# Missing values are handled without dropping samples:
#   Categorical NAs  -> modeled as an explicit "Missing" level
#   Continuous NAs   -> median-imputed, with a companion
#                       (1|<var>_missing) random effect that
#                       captures any systematic expression
#                       difference between observed and
#                       missing samples
# ============================================================

# ------------------------------------------------------------
# Classify each column as continuous or categorical.
# Numeric columns with few distinct integer values (e.g. Batch
# coded 1-8) are treated as categorical so they are not
# mistakenly modeled as dose-response variables.
# ------------------------------------------------------------
identify_column_types <- function(metadata, max_unique_for_categorical = 10) {
  cont_vars <- character()
  cat_vars  <- character()
  for (v in names(metadata)) {
    x     <- metadata[[v]]
    x_obs <- if (is.numeric(x)) x[!is.na(x)] else x
    looks_categorical <- is.numeric(x) &&
      length(x_obs) > 0 &&
      n_distinct(x_obs) <= max_unique_for_categorical &&
      all(x_obs == round(x_obs))
    if (!is.numeric(x) || looks_categorical)
      cat_vars <- c(cat_vars, v)
    else
      cont_vars <- c(cont_vars, v)
  }
  list(cont_vars = cont_vars, cat_vars = cat_vars)
}

# ------------------------------------------------------------
# Prepare metadata for modeling.
#
# Categorical NAs become an explicit "Missing" factor level so
# all samples remain in the model and the missingness pattern
# is captured by the random effect for that variable.
#
# Continuous NAs are median-imputed, and a binary indicator
# column (<var>_missing) is added as a random effect to absorb
# any expression difference between observed and imputed samples.
# ------------------------------------------------------------
prepare_metadata <- function(metadata, max_unique_for_categorical = 10) {
  types     <- identify_column_types(metadata, max_unique_for_categorical)
  cont_vars <- types$cont_vars
  cat_vars  <- types$cat_vars
  
  metadata <- metadata %>%
    mutate(across(all_of(cat_vars), as.factor)) %>%
    mutate(across(all_of(cat_vars), ~ fct_na_value_to_level(.x, "Missing")))
  
  missing_indicator_vars <- character()
  for (v in cont_vars) {
    if (any(is.na(metadata[[v]]))) {
      ind_name          <- paste0(v, "_missing")
      metadata[[ind_name]] <- factor(
        is.na(metadata[[v]]),
        levels = c(FALSE, TRUE),
        labels = c("Observed", "Missing")
      )
      missing_indicator_vars <- c(missing_indicator_vars, ind_name)
      metadata[[v]][is.na(metadata[[v]])] <- median(metadata[[v]], na.rm = TRUE)
    }
  }
  
  list(
    metadata               = metadata,
    cont_vars              = cont_vars,
    cat_vars               = cat_vars,
    missing_indicator_vars = missing_indicator_vars
  )
}

# ------------------------------------------------------------
# Assign each column a role: fixed, random, or excluded.
# All categoricals are random effects because variancePartition
# requires either all-fixed or all-random for categoricals, and
# random effects are the appropriate choice here.
# ------------------------------------------------------------
classify_variables <- function(metadata, cont_vars, cat_vars,
                               id_like_frac = 0.9) {
  n    <- nrow(metadata)
  rows <- list()
  
  for (v in cont_vars) {
    x <- metadata[[v]]
    if (length(x) == 0 || all(is.na(x)) || sd(x, na.rm = TRUE) == 0)
      rows[[v]] <- tibble(variable = v, role = "excluded", reason = "zero variance")
    else
      rows[[v]] <- tibble(variable = v, role = "fixed",    reason = "continuous")
  }
  
  for (v in cat_vars) {
    x <- droplevels(as.factor(metadata[[v]]))
    x <- x[x != "Missing"]
    k <- nlevels(x)
    if (k <= 1)
      rows[[v]] <- tibble(variable = v, role = "excluded", reason = "zero variance")
    else if (k / n > id_like_frac)
      rows[[v]] <- tibble(variable = v, role = "excluded", reason = "identifier-like")
    else
      rows[[v]] <- tibble(variable = v, role = "random",   reason = sprintf("%d levels", k))
  }
  
  bind_rows(rows)
}

# ------------------------------------------------------------
# Build the model formula from the variable classification.
# ------------------------------------------------------------
build_formula <- function(classification) {
  fixed_terms  <- classification$variable[classification$role == "fixed"]
  random_terms <- sprintf("(1 | %s)",
                          classification$variable[classification$role == "random"])
  rhs <- paste(c(fixed_terms, random_terms), collapse = " + ")
  if (nchar(rhs) == 0) return(NULL)
  as.formula(paste("~", rhs))
}

# Drop continuous fixed-effect variables that would cause rank
# deficiency in the design matrix. This can happen when a variable
# has zero or near-zero variance after imputation, or is an exact
# linear combination of other variables. Dropped variables are
# added to the excluded list with an explanatory reason.
drop_collinear_fixed <- function(expr_mat, classification, meta) {
  fixed_vars <- classification$variable[classification$role == "fixed"]
  
  if (length(fixed_vars) < 2) return(classification)
  
  X <- model.matrix(as.formula(paste("~", paste(fixed_vars, collapse = " + "))),
                    data = meta)
  X <- X[, -1, drop = FALSE]  # drop intercept
  
  rank_X    <- qr(X)$rank
  n_cols    <- ncol(X)
  
  if (rank_X == n_cols) return(classification)  # full rank, nothing to do
  
  # Identify redundant columns via pivoted QR decomposition and
  # flag the corresponding variables as excluded.
  pivot     <- qr(X)$pivot
  keep_cols <- pivot[seq_len(rank_X)]
  keep_vars <- colnames(X)[keep_cols]
  
  # Map design matrix column names back to variable names
  # (model.matrix may expand a variable name with a suffix).
  redundant_vars <- fixed_vars[!sapply(fixed_vars, function(v) any(startsWith(keep_vars, v)))]
  
  classification$role[classification$variable %in% redundant_vars]   <- "excluded"
  classification$reason[classification$variable %in% redundant_vars] <- "collinear with another fixed-effect variable"
  
  classification
}

round_number <- function(x) {
  ifelse(
    abs(x) < 0.001,
    format(signif(x, 2), scientific = TRUE),
    format(round(x, 3), nsmall = 3, trim = TRUE)
  )
}

# Returns a tibble of pairwise canonical correlations between all
# modeled variables, without duplicates (A vs B appears once only).
# Variables ending in "_missing" are excluded from CCA since their
# correlation with other variables reflects imputation structure
# rather than a real relationship between measured quantities.
# Each pair is computed using only samples where both variables are
# non-missing in the original data. Pairs with fewer than
# min_complete such samples receive NA.
compute_canonical_correlations <- function(form, meta, original_meta,
                                            min_complete = 10) {
  modeled_vars     <- all.vars(form)
  non_missing_vars <- modeled_vars[!endsWith(modeled_vars, "_missing")]

  if (length(non_missing_vars) < 2) {
    return(tibble(
      variable_1     = character(),
      variable_2     = character(),
      canonical_cor  = numeric()
    ))
  }

  get_mm <- function(v, data) {
    mm <- model.matrix(as.formula(paste("~ 0 +", v)), data = data)
    mm[, apply(mm, 2, var) > 0, drop = FALSE]
  }

  rows <- list()
  n_vars <- length(non_missing_vars)

  for (i in seq_len(n_vars - 1)) {
    for (j in seq(i + 1, n_vars)) {
      v1 <- non_missing_vars[i]
      v2 <- non_missing_vars[j]

      orig_cols <- intersect(c(v1, v2), names(original_meta))
      if (length(orig_cols) < 2) next
      complete <- complete.cases(original_meta[, orig_cols, drop = FALSE])

      cc <- if (sum(complete) < min_complete) {
        NA_real_
      } else {
        X1 <- tryCatch(get_mm(v1, meta[complete, , drop = FALSE]), error = function(e) NULL)
        X2 <- tryCatch(get_mm(v2, meta[complete, , drop = FALSE]), error = function(e) NULL)
        if (is.null(X1) || is.null(X2) || ncol(X1) == 0 || ncol(X2) == 0) {
          NA_real_
        } else {
          tryCatch(round(cancor(X1, X2)$cor[1], 3), error = function(e) NA_real_)
        }
      }

      rows[[length(rows) + 1]] <- tibble(
        variable_1    = v1,
        variable_2    = v2,
        canonical_cor = cc
      )
    }
  }

  bind_rows(rows)
}

# Exclude the alphabetically second variable from each pair with CCA = 1.
exclude_perfect_cca_variables <- function(classification, cca_tbl) {
  perfect <- cca_tbl %>%
    filter(!is.na(canonical_cor), canonical_cor == 1)

  if (nrow(perfect) == 0) {
    return(classification)
  }

  drop_info <- perfect %>%
    rowwise() %>%
    mutate(
      kept_var = sort(c(variable_1, variable_2))[1],
      drop_var = sort(c(variable_1, variable_2))[2]
    ) %>%
    ungroup() %>%
    distinct(drop_var, kept_var)

  for (i in seq_len(nrow(drop_info))) {
    v <- drop_info$drop_var[i]
    partner <- drop_info$kept_var[i]
    idx <- which(classification$variable == v)
    if (length(idx) == 1 && classification$role[idx] != "excluded") {
      classification$role[idx] <- "excluded"
      classification$reason[idx] <- paste0(
        "perfect canonical correlation (CCA = 1) with ", partner
      )
    }
  }

  classification
}

# ------------------------------------------------------------
# Run variance partitioning for one dataset.
#
# Returns a tibble with one row per metadata variable (including
# any <var>_missing indicators), plus an "Unexplained variance"
# row. Excluded variables have variance_explained = 0.
# ------------------------------------------------------------
run_variance_partition <- function(expr_mat, metadata,
                                   id_like_frac = 0.9,
                                   max_unique_for_categorical = 10,
                                   min_samples = 6,
                                   n_genes = 1000,
                                   random_seed = 1) {
  stopifnot(ncol(expr_mat) == nrow(metadata))
  
  if (ncol(expr_mat) < min_samples)
    stop(sprintf(
      "Only %d samples (min_samples = %d). Increase min_samples or check your data.",
      ncol(expr_mat), min_samples
    ))
  
  prepped <- prepare_metadata(metadata, max_unique_for_categorical)
  meta    <- prepped$metadata
  
  cls <- classify_variables(
    meta,
    cont_vars = prepped$cont_vars,
    cat_vars  = c(prepped$cat_vars, prepped$missing_indicator_vars),
    id_like_frac = id_like_frac
  )

  form_for_cca <- build_formula(cls)
  cor_tbl <- if (is.null(form_for_cca)) {
    tibble(
      variable_1    = character(),
      variable_2    = character(),
      canonical_cor = numeric()
    )
  } else {
    compute_canonical_correlations(form_for_cca, meta, original_meta = metadata) %>%
      arrange(desc(canonical_cor))
  }

  cls <- exclude_perfect_cca_variables(cls, cor_tbl)
  cls <- drop_collinear_fixed(expr_mat, cls, meta)

  excluded_tbl <- cls %>%
    filter(role == "excluded") %>%
    transmute(variable, role,
              variance_explained     = 0,
              status = paste0("excluded: ", reason))

  form <- build_formula(cls)
  if (is.null(form)) {
    return(list(
      variance_explained = excluded_tbl %>%
        arrange(desc(variance_explained)) %>%
        mutate(variance_explained = round_number(variance_explained)),
      cca = cor_tbl %>%
        mutate(canonical_cor = round_number(canonical_cor))
    ))
  }

  # This prevents errors for some datasets where metadata variables are on very different scales.
  fixed_vars <- cls$variable[cls$role == "fixed"]
  meta[, fixed_vars] <- scale(meta[, fixed_vars, drop = FALSE]) 

  n_cores <- 40
  BPPARAM <- if (n_cores > 1) MulticoreParam(workers = n_cores) else SerialParam()
#  cl <- makeCluster(n_cores)
#  registerDoParallel(cl)

  gene_var <- apply(expr_mat, 1, stats::var)
  varying_genes <- names(gene_var)[!is.na(gene_var) & gene_var > 0]
  if (length(varying_genes) == 0) {
    stop("No genes with non-zero variance remain after filtering.")
  }

  set.seed(random_seed)
  selected_genes <- if (length(varying_genes) > n_genes) {
    sample(varying_genes, n_genes)
  } else {
    varying_genes
  }
  expr_subset <- expr_mat[selected_genes, , drop = FALSE]

  print(paste0(
    "Fitting variance partitioning model on ",
    nrow(expr_subset), " of ", nrow(expr_mat), " genes with ",
    n_cores, " cores"
  ))
  varPart <- fitExtractVarPartModel(expr_subset, form, meta, BPPARAM = BPPARAM)
  print(paste0("Done fitting variance partitioning model"))

#  stopCluster(cl)

  varPart <- as.data.frame(varPart)
  varPart$feature <- rownames(varPart)
  
  fitted_tbl <- varPart %>%
    pivot_longer(-feature, names_to = "variable", values_to = "variance_explained") %>%
    group_by(variable) %>%
    summarise(variance_explained = mean(variance_explained), .groups = "drop") %>%
    mutate(
      variable = if_else(variable == "Residuals", "Unexplained variance", variable),
      role = case_when(
        variable == "Unexplained variance"              ~ "unexplained",
        variable %in% cls$variable[cls$role == "random"] ~ "random",
        TRUE                                            ~ "fixed"
      ),
      status = "ok"
    )

  list(
    variance_explained = bind_rows(fitted_tbl, excluded_tbl) %>%
#      mutate(n_samples = ncol(expr_mat)) %>%
#      relocate(n_samples, .before = variable) %>%
      arrange(desc(variance_explained)) %>%
      mutate(variance_explained = round_number(variance_explained)),
    cca = cor_tbl %>%
      mutate(canonical_cor = round_number(canonical_cor))
  ) %>%
    return()
}

processDataset <- function(dataset_id, expr_file_path, metadata_file_path, is_microarray, out_variance_file_path, out_cca_file_path) {
  if (file.exists(out_cca_file_path)) {
    return(NULL)
  }

  expr_data <- getExprData(expr_file_path)

  metadata <- getMetadata(metadata_file_path)
  platform_id <- metadata$Platform_ID
  metadata <- metadata$Metadata

  # These variables do not need to be included in this.
  #   Removing them speeds this analysis up.
  if (dataset_id %in% c("GSE62944_Normal", "GSE62944_Tumor")) {
    metadata <- select(metadata, -starts_with("icd_"), -all_of(c("ajcc_staging_edition", "bcr_patient_uuid")))
  }

  if (is.null(metadata)) {
    file.create(c(), out_variance_file_path)
    file.create(c(), out_cca_file_path)
    return(NULL)
  }

  sample_ids <- sort(intersect(colnames(expr_data), rownames(metadata)))

  if (length(sample_ids) < 5) {
    stop(paste0("There are few, if any, matching samples between the metadata and expression data for ", file_path1, "."))
  }

  expr_data <- expr_data[,sample_ids]
  metadata <- metadata[sample_ids, , drop=FALSE]

  print(paste0("Partitioning variance for ", expr_file_path))
  result <- run_variance_partition(expr_data, metadata)

  write_tsv(result$variance_explained, out_variance_file_path)
  write_tsv(result$cca, out_cca_file_path)
}

# Run dopplegangR for pairwise comparisons of datasets.
expr_file_paths <- list.files(datadir, full.names = T)

#################
#expr_file_paths <- expr_file_paths[grepl("METABRIC", expr_file_paths)]
#expr_file_paths <- expr_file_paths[grepl("GSE62944_Normal", expr_file_paths)]
expr_file_paths <- expr_file_paths[grepl("GSE62944_Tumor", expr_file_paths)]
#expr_file_paths <- expr_file_paths[grepl("GSE2990", expr_file_paths)]

#sequencing_platforms <- c("GPL18573", "GPL9052", "GPL11154", "GPL1791")

tmps <- c()
for (i in 1:length(expr_file_paths)) {
  dataset_id <- sub(".tsv.gz", "", basename(expr_file_paths[i]))

  expr_file_path <- expr_file_paths[i]
  metadata_file_path <- str_c(metadata_dir, "/", dataset_id, ".tsv")
  out_variance_file_path <- str_c(out_dir, "/", dataset_id, "_variance.tsv.gz")
  out_cca_file_path <- str_c(out_dir, "/", dataset_id, "_cca.tsv.gz")

  processDataset(dataset_id, expr_file_path, metadata_file_path, is_microarray, out_variance_file_path, out_cca_file_path)
#break
}
