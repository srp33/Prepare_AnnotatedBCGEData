
SCAN_normalise <- function(gseID, annotation_package, probe_summary, celFilePaths) {
    
    out_file_path <- paste0(normalized_data, gseID, ".tsv.gz")
    tmp_dir <- paste0("/tmp/", gseID)

    if (file.exists(out_file_path)) {
        print(paste0(gseID, " has already been processed!"))
    } else {
        unlink(tmp_dir, recursive = TRUE, force = TRUE)
        dir.create(tmp_dir)

        print(paste0("Downloading ", gseID, " for processing!"))
        GSE <- getGEOSuppFiles(gseID, filter_regex =  paste0(gseID, "_RAW.tar"))
        tmp <- rownames(GSE)
        untar(tmp[1], exdir = tmp_dir)

        all_normalized <- NULL

        for (celFilePath in celFilePaths) {
            gsm_id <- basename(celFilePath)
            gsm_id <- gsub("\\.cel.gz", "", gsm_id, ignore.case = TRUE)

            normalized <- exprs(SCAN(celFilePath, annotationPackageName = annotation_package, probeSummaryPackage = probe_summary))
            normalized <- as_tibble(normalized, rownames = "Gene")

            colnames(normalized)[2] <- gsm_id

            if (is.null(all_normalized)) {
                all_normalized <- normalized
            } else {
                all_normalized <- inner_join(all_normalized, normalized, by = "Gene")
            }
        }
        write_tsv(all_normalized, out_file_path)
        print(paste0("Saved to ", out_file_path))
    }
}


GSE23720 <- getGEO("GSE23720")
GSE23720_celfile <- pData(GSE23720[[1]]) %>%
  mutate(geo_accession = str_c("/tmp/GSE23720/", geo_accession, ".CEL.gz")) %>%
  pull(2)

#format to run function
#SCAN_normalise <- function(gseID, annotation_package, probe_summary, celFilePaths)

SCAN_normalise("GSE23720", "pd.hg.u133.plus.2", "hgu133plus2hsentrezgprobe", GSE23720_celfile)

unlink("GSE23720", recursive = TRUE, force = TRUE)
