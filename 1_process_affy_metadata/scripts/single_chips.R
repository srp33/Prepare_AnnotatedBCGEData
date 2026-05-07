#This script parses gene expression metadata

gseIDs <- read_tsv("/Data/gseIDs.tsv", comment = "#")


keep_varible <- c("Affymetrix Human Exon 1.0 ST Array [transcript (gene) version]", "Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]",
                  "Affymetrix Human Genome U95 Version 2 Array", "Affymetrix GeneChip HT-HG_U133A Early Access Array",
                  "Affymetrix Human Genome U133A Array", "Affymetrix Human Genome U133A 2.0 Array", "Affymetrix Human Genome U133 Plus 2.0 Array")


gseID_list <- NULL

#create a new list with gseID's we want to keep
for (i in seq_along(gseIDs$geneChip)) {
  keep_sample <- gseIDs$geneChip[i]
  if (keep_sample %in% keep_varible) {
    gseID_list <- rbind(gseID_list, gseIDs$gseID[i])
  }
}

for (gseID in gseID_list) {

  out_file_path_metadata <- paste0(metadata_dir, gseID, ".tsv")
  out_file_path_raw_metadata <- paste0(raw_metadata_dir, gseID, ".tsv")

  if (file.exists(out_file_path_metadata)) {
      print(paste0(gseID, " has already been processed!"))
  } else {
      # write un-curated metadata to file
      df <- getFromGEO(gseID)      
      write_tsv(df, file.path(raw_metadata_dir, paste0(gseID, ".tsv")))

      #some housekeeping and cleaning up column and column names
      metadata <- df %>%
        removeUnusefulCols() %>%
        rename_with(~str_replace_all(., "_ch1", "")) %>%
        mutate(across(where(is.character), ~replace(., . %in% c("?", "--", ""), NA)))

      # source file to address special cases
      source("scripts/special_cases.R")
      
      # summarise variables
      varSummary <- summariseVariables(metadata)

      #write cleaned up data to files
      writeOutput(gseID)
  }
}
