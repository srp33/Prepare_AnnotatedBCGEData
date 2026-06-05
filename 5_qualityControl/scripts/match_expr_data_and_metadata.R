expr_dir <- "/Data/expression_data2"
file_paths_expr <- list.files(expr_dir, full.names = T)

hgnc_special_cases <- c("GSE62944_Tumor", "GSE62944_Normal", "GSE81538", "GSE96058_HiSeq", "GSE96058_NextSeq", "ICGC_KR", "METABRIC")

ensg_special_cases <- c("SCANB.9206", "OSLO2EMIT0.103", "Normal.66", "ABiM.405", "ABiM.100")

for (expr_file_path in file_paths_expr) {
    dataset_id <- basename(expr_file_path)
    dataset_id <- sub(".tsv.gz", "", dataset_id)

    meta_file_path = paste0("/Data/prelim_metadata/", dataset_id, ".tsv")

    expr_data <- read_tsv(expr_file_path)
    metadata <- read_tsv(meta_file_path)

    out_expr_path <- paste0("/Data/expression_data3/", basename(expr_file_path))
    out_meta_path <- paste0("/Data/prelim_metadata2/", basename(meta_file_path))

    if (file.exists(out_expr_path) & file.exists(out_meta_path)) {
      print(paste0(out_expr_path, " and ", out_meta_path, " already exist."))
      next
    }

    Sample_ID <- names(expr_data)[3:ncol(expr_data)]
    all_samples <- metadata %>%
        pull(Sample_ID)
    keep_samples <- sort(intersect(Sample_ID, all_samples))

    clean_metadata <- arrange(metadata, match(Sample_ID, keep_samples))

    if (dataset_id %in% hgnc_special_cases) {
      clean_expr_data <- dplyr::select(expr_data, Dataset_ID, HGNC_Symbol, all_of(keep_samples))
    } else {
        if (dataset_id %in% ensg_special_cases) {
          clean_expr_data <- dplyr::select(expr_data, Dataset_ID, ENSG, all_of(keep_samples))
        } else {
          clean_expr_data <- dplyr::select(expr_data, Dataset_ID, Gene, all_of(keep_samples))
        }
    }

    print(paste0("Writing to ", out_expr_path, " and ", out_meta_path, "!"))

    write_tsv(clean_expr_data, out_expr_path)
    write_tsv(clean_metadata, out_meta_path)
}
