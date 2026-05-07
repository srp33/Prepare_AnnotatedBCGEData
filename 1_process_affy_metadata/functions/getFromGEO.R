# This script downloads data from GEO, selects the metadata, and cleans up the medadata.

getFromGEO <- function(geoID) {
#   geoID <- "GSE22093"
  gseData <- getGEO(geoID)
  df <- gseData[[1]]
  
  meta <- pData(df) |>
    clean_names() |>
    dplyr::rename(Sample_ID = geo_accession) |>
    dplyr::select(Sample_ID, everything()) |> 
    mutate(Dataset_ID = geoID, .before = Sample_ID) |>
    mutate(Platform_ID = platform_id, .after = Sample_ID) |>
    dplyr::select(-platform_id) 
    
  meta <- replace(meta, meta=='n/a', NA)
  meta <- replace(meta, meta=='N/A', NA)

  return(metadata = meta)
}