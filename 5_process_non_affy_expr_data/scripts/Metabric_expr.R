
# Cbioportal hub
# https://github.com/cBioPortal/datahub/tree/master/public/brca_metabric

# download the patient data from this URL
download.file("https://osf.io/download/x8hbc/", 
              destfile = paste0(tmp_dir, "metabric_expr.txt"))

metaBric <- read_tsv(paste0(tmp_dir, "metabric_expr.txt"), comment = "#") %>%
  rename("HGNC_Symbol" = "Hugo_Symbol") %>%
  dplyr::select(-Entrez_Gene_Id)

print("Writing METABRIC to file!")
write_tsv(metaBric, paste0(data_dir, "METABRIC.tsv.gz"))



