file_paths <- list.files("/Data/expression_data3", full.names = T)

hgnc_special_cases <- c("GSE62944_Tumor.tsv.gz", "GSE62944_Normal.tsv.gz", "GSE81538.tsv.gz", "GSE96058_HiSeq.tsv.gz",
                        "GSE96058_NextSeq.tsv.gz", "ICGC_KR.tsv.gz", "METABRIC.tsv.gz")

ensg_special_cases <- c("SCANB.9206.tsv.gz", "OSLO2EMIT0.103.tsv.gz", "Normal.66.tsv.gz", "ABiM.100.tsv.gz", "ABiM.405.tsv.gz")

# Retrieve and fine tune gene annotations
download.file("https://zenodo.org/records/20549560/files/HGNC_mappings.tsv.gz?download=1", "/tmp/gene_mappings.tsv.gz")
gene_annotations <- read_tsv("/tmp/gene_mappings.tsv.gz", comment = "#") %>%
  dplyr::rename(Ensembl_Gene_ID = `Ensembl gene ID`,
         Entrez_Gene_ID = `NCBI gene ID`,
         HGNC_Symbol = `Approved symbol`,
         Chromosomal_Band = `Chromosome location`) %>%
    filter(`Locus group` == "protein-coding gene") %>%
    filter(!is.na(Ensembl_Gene_ID)) %>%
    filter(!is.na(Entrez_Gene_ID)) %>%
    filter(!is.na(HGNC_Symbol)) %>%
    mutate(Entrez_Gene_ID = as.integer(Entrez_Gene_ID)) %>%
    dplyr::select(Ensembl_Gene_ID, Entrez_Gene_ID, HGNC_Symbol, Chromosomal_Band)

for (file_path in file_paths) {
  cat("\n")
  print(paste0("Reading in ", file_path, "!"))
  cat("\n")

  file_name <- file_path %>% basename()
  out_file_path <- paste0("/Data/expression_data4/", file_name)

  if (file.exists(out_file_path)) {
    print(paste0(out_file_path, " already exists."))
    next
  }

  expr_data <- read_tsv(file_path)

  if (file_name %in% hgnc_special_cases) {
      join_column <- "HGNC_Symbol"
  } else {
    if (file_name %in% ensg_special_cases) {
      expr_data <- dplyr::rename(expr_data, Ensembl_Gene_ID = ENSG) %>%
        filter(str_detect(Ensembl_Gene_ID, "_PAR_Y$", negate = TRUE)) %>%
        mutate(Ensembl_Gene_ID = str_extract(Ensembl_Gene_ID, "ENSG\\d+"))

      join_column <- "Ensembl_Gene_ID"
    } else {
        expr_data <- mutate(expr_data, Gene = str_replace(Gene, "_at", "")) %>%
          dplyr::rename(Entrez_Gene_ID = Gene) %>%
          filter(!str_detect(Entrez_Gene_ID, "^AFFX")) %>%
          mutate(Entrez_Gene_ID = as.integer(Entrez_Gene_ID))

        join_column <- "Entrez_Gene_ID"
    }
  }

  expr_data <- inner_join(gene_annotations, expr_data, by = join_column) %>%
    dplyr::select(Dataset_ID, Entrez_Gene_ID, Ensembl_Gene_ID, HGNC_Symbol, Chromosomal_Band, everything())

  group_by(expr_data, Ensembl_Gene_ID) %>% dplyr::summarize(Count = n()) %>% dplyr::filter(Count > 1) %>% print()
  group_by(expr_data, Entrez_Gene_ID) %>% dplyr::summarize(Count = n()) %>% dplyr::filter(Count > 1) %>% print()
  group_by(expr_data, HGNC_Symbol) %>% dplyr::summarize(Count = n()) %>% dplyr::filter(Count > 1) %>% print()

  expr_data <- arrange(expr_data, Entrez_Gene_ID)
 
  write_tsv(expr_data, out_file_path)
}
