
SCAN_normalise <- function(gseID, annotation_package, probe_summary, GSM_to_exclude = c()) {

    out_file_path <- paste0(normalized_data, gseID, ".tsv.gz")
    tmp_dir <- paste0("/tmp/", gseID)
   
    if (file.exists(out_file_path)) {
        print(paste0(gseID, " has already been processed!"))
    } else {
        unlink(tmp_dir, recursive = TRUE, force = TRUE)
        dir.create(tmp_dir)

        print(paste0("Downloading ", gseID, " for processing!"))
        GSE <- getGEOSuppFiles(gseID, makeDirectory = F, baseDir = tmp_dir, filter_regex = "*.tar")
        tmp <- rownames(GSE)
        untar(tmp[1], exdir = tmp_dir)

        celFilePaths <- list.files(tmp_dir, pattern = "*.CEL", full.names = T, ignore.case = T)

        all_normalized <- NULL

        for (celFilePath in celFilePaths) {
            gsm_id <- basename(celFilePath)
            gsm_id <- gsub("\\.cel.gz", "", gsm_id, ignore.case = TRUE)

            if (gsm_id %in% GSM_to_exclude) {
                next
            }

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

        unlink(tmp_dir, recursive = TRUE, force = TRUE)
    }
}

# To run individual GSE ids, run the code line below, substituting with actual values (GSEID, annotation_package, probe_summary)
# e.g SCAN_normalise("GSE118432", "pd.hugene.1.0.st.v1", "hugene10sthsentrezgprobe")

for (gseID in huExon$gseID) {
SCAN_normalise(gseID, "pd.huex.1.0.st.v2", "huex10sthsentrezgprobe")
}

for (gseID in huGene$gseID) {
  SCAN_normalise(gseID, "pd.hugene.1.0.st.v1", "hugene10sthsentrezgprobe")
}

for (gseID in U95_2$gseID) {
  SCAN_normalise(gseID, "pd.hg.u95av2", "hgu95av2hsentrezgprobe")
}

for (gseID in U133_A_Early_Access$gseID) {
  SCAN_normalise(gseID, "pd.ht.hg.u133a", "u133aaofav2hsentrezgprobe")
}

for (gseID in U133_A$gseID) {
SCAN_normalise(gseID, "pd.hg.u133a", "hgu133ahsentrezgprobe")
}

for (gseID in U133_A2$gseID) {
SCAN_normalise(gseID, "pd.hg.u133a.2", "hgu133a2hsentrezgprobe")
}

# we are excluding c("GSM125119", "GSM125120") in GSE5460 because the are listed as corrupted on GEO website
for (gseID in U133_plus_2$gseID) {
    SCAN_normalise(gseID, "pd.hg.u133.plus.2", "hgu133plus2hsentrezgprobe", c("GSM125119", "GSM125120"))
}
