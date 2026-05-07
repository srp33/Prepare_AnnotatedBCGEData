library(GEOquery)
library(tidyverse)

dual_chips <- c("GSE1456", "GSE3494", "GSE4922")

getDualGEO <- function(geoID) {
  gseData <- getGEO(geoID)

  U133A <- gseData[[1]]
  U133B <- gseData[[2]]

  meta_A <- pData(U133A) %>%
    as_tibble(rownames = NA) %>%
    clean_names() %>%
    dplyr::rename(Sample_ID = geo_accession) %>%
    dplyr::select(Sample_ID, everything()) %>%
    mutate(Dataset_ID = (paste0(geoID, "_U133A")), .before = Sample_ID) %>%
    mutate(Platform_ID = platform_id, .after = Sample_ID) %>%
    dplyr::select(-platform_id)

  meta_B <- pData(U133B) %>%
    as_tibble(rownames = NA) %>%
    clean_names() %>%
    dplyr::rename(Sample_ID = geo_accession) %>%
    dplyr::select(Sample_ID, everything()) %>%
    mutate(Dataset_ID = (paste0(geoID, "_U133B")), .before = Sample_ID) %>%
    mutate(Platform_ID = platform_id, .after = Sample_ID) %>%
    dplyr::select(-platform_id)

  return(list(metadata_A = meta_A, metadata_B = meta_B))
}

clean_metadata_Dual <- function(meta) {
  metadata <- meta %>%
    removeUnusefulCols() %>%
    rename_with(~str_replace_all(., "_ch1", "")) %>%
    dplyr::select(-c("title", "description")) %>%
    dplyr::select(-starts_with("characteristics"))

    # remove colums with more than 50% NA
    # metadata <- metadata[, colMeans(is.na(metadata)) < 0.5]

  return(metadata)
}

for (gseID in dual_chips) {
    df <- getDualGEO(gseID)

    # write un-curated metadata to file
    write_tsv(df$metadata_A, file.path(raw_metadata_dir, paste0(gseID, "_U133A.tsv")))
    write_tsv(df$metadata_B, file.path(raw_metadata_dir, paste0(gseID, "_U133B.tsv")))

    # remove unuseful columns
    metadata_A <- clean_metadata_Dual(df$metadata_A)
    metadata_B <- clean_metadata_Dual(df$metadata_B)

    if (gseID == "GSE4922") {
      metadata_A <- metadata_A %>%
        dplyr::select(-c("all_patients_1_included_in_survival_analysis",
                         "er_endocrine_therapy_only_1_included_in_survival_analysis", 
                         "genetic_grade_signature_status_prediction_by_sws_classifier",
                         "no_systemic_therapy_1_included_in_survival_analysis",
                         "probability_1_like_by_sws_classifier",
                         "probability_3_like_by_sws_classifier")) 

      metadata_B <- metadata_B %>%
        dplyr::select(-c("all_patients_1_included_in_survival_analysis",
                         "er_endocrine_therapy_only_1_included_in_survival_analysis", 
                         "genetic_grade_signature_status_prediction_by_sws_classifier",
                         "no_systemic_therapy_1_included_in_survival_analysis",
                         "probability_1_like_by_sws_classifier",
                         "probability_3_like_by_sws_classifier"))
    }

    # summarise variables
    varSummary_A <- summariseVariables(metadata_A)
    varSummary_B <- summariseVariables(metadata_B)

    #write cleaned up data to files
    if ((ncol(metadata_A) > 2)) {
      write_tsv(metadata_A, file.path(metadata_dir, paste0(gseID, "_U133A.tsv")))
      write_tsv(metadata_B, file.path(metadata_dir, paste0(gseID, "_U133B.tsv")))
    }

    if (!is.null(varSummary_A$numSummary)) {
      write_tsv(varSummary_A$numSummary, file.path(metadata_summaries, paste0(gseID, "_U133A_num.tsv")))
      write_tsv(varSummary_B$numSummary, file.path(metadata_summaries, paste0(gseID, "_U133B_num.tsv")))
    }

    if (!is.null(varSummary_A$charSummary)) {
      write_tsv(varSummary_A$charSummary, file.path(metadata_summaries, paste0(gseID, "_U133A_char.tsv")))
      write_tsv(varSummary_B$charSummary, file.path(metadata_summaries, paste0(gseID, "_U133B_char.tsv")))
    }
}

gseData <- getGEO("GSE6532")

U133A <- gseData[[2]]
U133B <- gseData[[3]]
U133Plus2 <- gseData[[1]]

parse_metadata <- function(chipID, platform) {
  meta <- pData(chipID) %>%
    clean_names() %>%
    dplyr::rename(Sample_ID = geo_accession) %>%
    dplyr::select(Sample_ID, everything()) %>%
    mutate(Dataset_ID = (paste0("GSE6532", platform)), .before = Sample_ID)%>%
    mutate(Platform_ID = platform_id, .after = Sample_ID) %>%
    dplyr::select(-platform_id)

  return(metadata = meta)
}

df_U133A <- parse_metadata(U133A, "_U133A")
df_U133B <- parse_metadata(U133B, "_U133B")
df_U133Plus2 <- parse_metadata(U133Plus2, "_U133Plus2")

# write un-curated metadata to file
write_tsv(df_U133A, file.path(raw_metadata_dir, paste0("GSE6532_U133A.tsv")))
write_tsv(df_U133B, file.path(raw_metadata_dir, paste0("GSE6532_U133B.tsv")))
write_tsv(df_U133Plus2, file.path(raw_metadata_dir, paste0("GSE6532_U133Plus2.tsv")))

clean_metadata <- function(chip) {
  metadata <- chip %>%
    removeUnusefulCols() %>%
    rename_with(~str_replace_all(., "_ch1", "")) %>%
    dplyr::select(-c("title", "description")) %>%
    mutate(across(where(is.character), ~replace(., . %in% c("KJ67", "KJ68", "KJ69", "KJX46", "KJX38", "KJ117"), NA)))

    # remove colums with more than 50% NA
    # metadata <- metadata[, colMeans(is.na(metadata)) < 0.5]

  return(metadata)
}

metadata_A <- clean_metadata(df_U133A) %>%      #The following variables are present in A but not B & C (distant_rfs, ggi, time_rfs)
  dplyr::select(-starts_with(c("characteristics_", "ggi")))
metadata_B <- clean_metadata(df_U133B)
metadata_C <- clean_metadata(df_U133Plus2)

# summarise variables
varSummary_A <- summariseVariables(metadata_A)
varSummary_B <- summariseVariables(metadata_B)
varSummary_C <- summariseVariables(metadata_C)

#write cleaned up data to files
write_tsv(metadata_A, file.path(metadata_dir, paste0("GSE6532_U133A.tsv")))
write_tsv(metadata_B, file.path(metadata_dir, paste0("GSE6532_U133B.tsv")))
write_tsv(metadata_C, file.path(metadata_dir, paste0("GSE6532_U133Plus2.tsv")))

if (!is.null(varSummary_A$numSummary)) {
  write_tsv(varSummary_A$numSummary, file.path(metadata_summaries, paste0("GSE6532_U133A_num.tsv")))
  write_tsv(varSummary_B$numSummary, file.path(metadata_summaries, paste0("GSE6532_U133B_num.tsv")))
  write_tsv(varSummary_C$numSummary, file.path(metadata_summaries, paste0("GSE6532_U133Plus2_num.tsv")))
}

if (!is.null(varSummary_A$charSummary)) {
  write_tsv(varSummary_A$charSummary, file.path(metadata_summaries, paste0("GSE6532_U133A_char.tsv")))
  write_tsv(varSummary_B$charSummary, file.path(metadata_summaries, paste0("GSE6532_U133B_char.tsv")))
  write_tsv(varSummary_C$charSummary, file.path(metadata_summaries, paste0("GSE6532_U133Plus2_char.tsv")))
}