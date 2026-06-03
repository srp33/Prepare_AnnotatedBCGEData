# Cbioportal hub
# https://github.com/cBioPortal/datahub/tree/master/public/brca_metabric

# download the patient data from this URL
download.file("https://zenodo.org/records/20097812/files/METABRIC_data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt.gz?download=1", 
              destfile = paste0(tmp_dir, "metabric_expr.txt"))

metaBric <- read_tsv(paste0(tmp_dir, "metabric_expr.txt"), comment = "#") %>%
  rename("HGNC_Symbol" = "Hugo_Symbol") %>%
  dplyr::select(-Entrez_Gene_Id)

print("Writing METABRIC to file!")
write_tsv(metaBric, paste0(data_dir, "METABRIC.tsv.gz"))
