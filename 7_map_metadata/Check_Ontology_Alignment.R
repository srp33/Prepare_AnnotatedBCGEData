library(tidyverse)
library(readxl)

# Checking against version 26.04d of NCIT
tmp_file_path = paste0(tempdir(), "/NCIT.csv.gz")
if (!file.exists(tmp_file_path)) {
  download.file("https://data.bioontology.org/ontologies/NCIT/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb&download_format=csv", tmp_file_path)
}
ncit = read_csv(tmp_file_path) %>%
  dplyr::select(`Preferred Label`, code) %>%
  dplyr::rename(Preferred = `Preferred Label`) %>%
  mutate(Together = str_c(Preferred, " (Code ", code, ")")) %>%
  pull(Together)

numeric = read_xlsx("Metadata_Mappings_Deduplicated.xlsx", sheet = "Numeric", col_types = rep("text", 7)) %>%
  separate_longer_delim(NCIT_field, delim = "||") %>%
  separate_longer_delim(NCIT_values, delim = "||")

filter(numeric, !(NCIT_field %in% ncit)) %>%
  View()

filter(numeric, !(NCIT_values %in% ncit)) %>%
  View()

categorical = read_xlsx("Metadata_Mappings_Deduplicated.xlsx", sheet = "Categorical", col_types = rep("text", 7)) %>%
  separate_longer_delim(NCIT_field, delim = "||") %>%
  separate_longer_delim(NCIT_values, delim = "||") %>%
  mutate(NCIT_field = str_trim(NCIT_field)) %>%
  mutate(NCIT_values = str_trim(NCIT_values))

filter(categorical, !(NCIT_field %in% ncit)) %>%
  View()

# filter(categorical, !(NCIT_field %in% ncit)) %>%
#   pull(NCIT_field) %>%
#   unique() %>%
#   sort() %>%
#   print()

filter(categorical, !(NCIT_values %in% ncit)) %>%
  View()

# filter(categorical, !(NCIT_values %in% ncit)) %>%
#   pull(NCIT_values) %>%
#   unique() %>%
#   sort() %>%
#   print()