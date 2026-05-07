
SCAN_normalise <- function(gseID, annotation_package, probe_summary, output_filename, celFilePaths) {
    out_file_path <- paste0(normalized_data, output_filename, ".tsv.gz")
    tmp_dir <- paste0("/tmp/", output_filename)

    if (file.exists(out_file_path)) {
        print(paste0(output_filename, " has already been processed!"))
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

#process GSE1456
GSE1456 <- getGEO("GSE1456")
GSE1456_U133A_celfile <- pData(GSE1456[[1]]) %>%
  mutate(geo_accession = str_c("/tmp/GSE1456_U133A/", geo_accession, ".CEL.gz")) %>%
  pull(2)
GSE1456_U133B_celfile  <- pData(GSE1456[[2]])  %>%
  mutate(geo_accession = str_c("/tmp/GSE1456_U133B/", geo_accession, ".CEL.gz")) %>%
  pull(2)

unlink("GSE1456", recursive = TRUE, force = TRUE)

#process GSE3494
GSE3494_U133A_celfile <- getGSEDataTables("GSE3494")[[2]] %>%
filter(`Affy platform` == "HG-U133A")  %>%
  mutate(`GEO Sample Accession #` = str_c("/tmp/GSE3494_U133A/", `GEO Sample Accession #`, ".CEL.gz")) %>%
  pull(1)

GSE3494_U133B_celfile <- getGSEDataTables("GSE3494")[[2]] %>%
  filter(`Affy platform` == "HG-U133B")  %>%
  mutate(`GEO Sample Accession #` = str_c("/tmp/GSE3494_U133B/", `GEO Sample Accession #`, ".CEL.gz")) %>%
  pull(1)

unlink("GSE3494", recursive = TRUE, force = TRUE)

#process GSE4922
GSE <- getGEOSuppFiles("GSE4922", filter_regex = "GSE4922_Clinical_file_for_both_Uppsala_Singapore_Samples.txt")
GSE4922 <- rownames(GSE) %>%
  read_tsv(col_names = TRUE) %>%
  separate(1, into = c("gsmID_U133A", "gsmID_U133B"), sep = "/") %>%
  mutate(gsmID_U133A = str_c("/tmp/GSE4922_U133A/", gsmID_U133A, ".CEL.gz"))  %>%
  mutate(gsmID_U133B = str_c("/tmp/GSE4922_U133B/", gsmID_U133B, ".CEL.gz"))

GSE4922_U133A_celfile <- pull(GSE4922, 1)
GSE4922_U133B_celfile <- pull(GSE4922, 2)

unlink("GSE4922", recursive = TRUE, force = TRUE)

#process GSE6532
GSE6532 <- getGEOSuppFiles("GSE6532", filter_regex = "GSE6532_LUMINAL_demo.txt.gz")
GSE6532_U133A_celfile <- rownames(GSE6532) %>%
  read_tsv(col_names = TRUE) %>%
  mutate(`geo_accn_hg-u133a` = str_c("/tmp/GSE6532_U133A/", `geo_accn_hg-u133a`, ".CEL.gz")) %>%
  pull(`geo_accn_hg-u133a`) %>%
  na.omit()

GSE6532_U133B_celfile <- rownames(GSE6532) %>%
  read_tsv(col_names = TRUE) %>%
  mutate(`geo_accn_hg-u133b` = str_c("/tmp/GSE6532_U133B/", `geo_accn_hg-u133b`, ".CEL.gz")) %>%
  pull(`geo_accn_hg-u133b`) %>%
  na.omit()

GSE6532_U133_2_celfile <- rownames(GSE6532) %>%
  read_tsv(col_names = TRUE) %>%
  mutate(`geo_accn_hg-u133plus2` = str_c("/tmp/GSE6532_U133Plus2/", `geo_accn_hg-u133plus2`, ".CEL.gz")) %>%
  pull(`geo_accn_hg-u133plus2`) %>%
  na.omit()

unlink("GSE6532", recursive = TRUE, force = TRUE)

#format to run function
#SCAN_normalise <- function(gseID, annotation_package, probe_summary, output_filename, celFilePaths)

SCAN_normalise("GSE1456", "pd.hg.u133a", "hgu133ahsentrezgprobe", "GSE1456_U133A", GSE1456_U133A_celfile)
SCAN_normalise("GSE1456", "pd.hg.u133b", "hgu133bhsentrezgprobe", "GSE1456_U133B", GSE1456_U133B_celfile)

SCAN_normalise("GSE3494", "pd.hg.u133a", "hgu133ahsentrezgprobe", "GSE3494_U133A", GSE3494_U133A_celfile)
SCAN_normalise("GSE3494", "pd.hg.u133b", "hgu133bhsentrezgprobe", "GSE3494_U133B", GSE3494_U133B_celfile)

SCAN_normalise("GSE4922", "pd.hg.u133a", "hgu133ahsentrezgprobe", "GSE4922_U133A", GSE4922_U133A_celfile)
SCAN_normalise("GSE4922", "pd.hg.u133b", "hgu133bhsentrezgprobe", "GSE4922_U133B", GSE4922_U133B_celfile)

SCAN_normalise("GSE6532", "pd.hg.u133a", "hgu133ahsentrezgprobe", "GSE6532_U133A", GSE6532_U133A_celfile)
SCAN_normalise("GSE6532", "pd.hg.u133b", "hgu133bhsentrezgprobe", "GSE6532_U133B", GSE6532_U133B_celfile)
SCAN_normalise("GSE6532", "pd.hg.u133.plus.2", "hgu133plus2hsentrezgprobe", "GSE6532_U133Plus2", GSE6532_U133_2_celfile)
