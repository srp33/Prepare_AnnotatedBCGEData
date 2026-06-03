library(tidyverse)
library(readxl)
library(writexl)

# Read from files and remove code prefixes and unnecessary columns.
numeric = read_xlsx("Metadata_Mappings_Deduplicated.xlsx", sheet = "Numeric", col_types = rep("text", 7)) %>%
  mutate(NCIT_field = str_trim(NCIT_field)) %>%
  mutate(NCIT_values = str_trim(NCIT_values)) %>%
  mutate(NCIT_field = str_replace_all(NCIT_field, " \\(Code ", " \\(")) %>%
  mutate(NCIT_values = str_replace_all(NCIT_values, " \\(Code ", " \\(")) %>%
  dplyr::select(-primitive_type)

categorical = read_xlsx("Metadata_Mappings_Deduplicated.xlsx", sheet = "Categorical", col_types = rep("text", 7)) %>%
  mutate(NCIT_field = str_trim(NCIT_field)) %>%
  mutate(NCIT_values = str_trim(NCIT_values)) %>%
  mutate(NCIT_field = str_replace_all(NCIT_field, " \\(Code ", " \\(")) %>%
  mutate(NCIT_values = str_replace_all(NCIT_values, " \\(Code ", " \\(")) %>%
  dplyr::rename(basic_type = primitive_type)

# Separate out date and identifier rows, add basic data types.
date = dplyr::filter(numeric, str_detect(NCIT_values, "Date \\(C25164\\)")) %>%
  mutate(basic_type = "Date Data Type (C48871)") %>%
  dplyr::select(-orig_values) %>%
  distinct()

numeric = dplyr::filter(numeric, str_detect(NCIT_values, "Date \\(C25164\\)", negate = TRUE)) %>%
  mutate(basic_type = "Numeric") %>%
  dplyr::select(-orig_values) %>%
  distinct()

identifier = dplyr::filter(categorical, str_detect(NCIT_values, "C25364")) %>%
  mutate(basic_type = "Instance Identifier Data Type (C95664)") %>%
  dplyr::select(-orig_values) %>%
  distinct()

categorical = dplyr::filter(categorical, str_detect(NCIT_values, "C25364", negate = TRUE)) %>%
  mutate(basic_type = "Discrete Set String Data Type (C95648)")

################################################################################

# Get the unique values for all fields.
all_fields = NULL

for (metadata_file_path in list.files(path = "../Data/prelim_metadata", full.names = TRUE)) {
  this_dataset = basename(metadata_file_path)
  this_dataset = sub(".tsv", "", this_dataset)
  
  metadata_dataset = read_tsv(metadata_file_path) %>%
    dplyr::select(-Dataset_ID, -Sample_ID, -Platform_ID)
  
  numeric_dataset = dplyr::filter(numeric, dataset == this_dataset)
  categorical_dataset = dplyr::filter(categorical, dataset == this_dataset)
  
  for (field in colnames(metadata_dataset)) {
    values = pull(metadata_dataset, all_of(field)) %>%
      unique() %>%
      sort(na.last = TRUE) %>%
      as.character()
    values[is.na(values)] <- "NA"
    
    fields = tibble(dataset = this_dataset, orig_field = field, orig_values = values)
    all_fields = bind_rows(all_fields, fields)
  }
}

# Remove the categorical fields because they are already separated out.
all_fields = anti_join(all_fields, categorical)

################################################################################

# Retrieve the unique original values.
identifier = inner_join(identifier, all_fields)
numeric = inner_join(numeric, all_fields) %>%
  mutate(orig_values = if_else(orig_values == "ND", "NA", orig_values)) %>%
  mutate(NCIT_values = if_else(dataset == "GSE20271" & orig_field == "post_chemo_size_cm" & orig_values == "Not medible", "Unknown (C17998)", NCIT_values)) %>%
  mutate(NCIT_values = if_else(dataset == "GSE20271" & orig_field == "post_chemo_size_cm" & orig_values == "residual cancer", "Unknown (C17998)", NCIT_values))
date = inner_join(date, all_fields) %>%
  mutate(orig_values = if_else(orig_values == "ND", "NA", orig_values))

# Combine all into one
# Re-order and rename columns.
# Re-sort.
# Find rows with NA in orig_values and update ontology term.
# Use Null basic data type for Unknown and NA values.
# Infer Integer vs. Float for all numeric values.
# Resave the Excel file under a different name.

out_file_path = "Metadata_Mappings_Final.xlsx"

all = bind_rows(numeric, date, identifier, categorical) %>%
  dplyr::select(dataset, orig_field, basic_type, NCIT_field, orig_values, NCIT_values, comments) %>%
  dplyr::rename(Dataset_ID = dataset,
                Original_Field_Name = orig_field,
                Field_Name_Ontology_Terms = NCIT_field,
                Original_Value = orig_values,
                Value_Ontology_Terms = NCIT_values,
                Basic_Data_Type_Ontology_Term = basic_type,
                Curatorial_Comments = comments) %>%
  arrange(Dataset_ID, Original_Field_Name, Basic_Data_Type_Ontology_Term, Field_Name_Ontology_Terms, Value_Ontology_Terms) %>%
  mutate(Dataset_ID = str_trim(Dataset_ID)) %>%
  mutate(Original_Field_Name = str_trim(Original_Field_Name)) %>%
  mutate(Field_Name_Ontology_Terms = str_trim(Field_Name_Ontology_Terms)) %>%
  mutate(Original_Value = str_trim(Original_Value)) %>%
  mutate(Value_Ontology_Terms = str_trim(Value_Ontology_Terms)) %>%
  mutate(Basic_Data_Type_Ontology_Term = str_trim(Basic_Data_Type_Ontology_Term)) %>%
  mutate(Curatorial_Comments = str_trim(Curatorial_Comments)) %>%
  mutate(Value_Ontology_Terms = if_else(Original_Value == "NA", "Not Applicable (C48660)", Value_Ontology_Terms)) %>%
  mutate(Basic_Data_Type_Ontology_Term = if_else(Value_Ontology_Terms == "Not Applicable (C48660)", "Null (C47840)", Basic_Data_Type_Ontology_Term)) %>%
  mutate(Basic_Data_Type_Ontology_Term = if_else(Value_Ontology_Terms == "Unknown (C17998)", "Null (C47840)", Basic_Data_Type_Ontology_Term)) %>%
  mutate(
    Basic_Data_Type_Ontology_Term = case_when(
      Basic_Data_Type_Ontology_Term == "Numeric" &
        str_detect(Original_Value, "^[+-]?\\d+$") ~
        "Integer Data Type (C95821)",
      
      Basic_Data_Type_Ontology_Term == "Numeric" &
        str_detect(Original_Value, "^[+-]?(\\d+\\.\\d*|\\.\\d+|\\d+)([eE][+-]?\\d+)?$") ~
        "Float (C48150)",

      Basic_Data_Type_Ontology_Term == "Numeric" &
        str_detect(Original_Value, "^\\d+/\\d+$") ~
        "Fraction (C25514)",
      
      TRUE ~ Basic_Data_Type_Ontology_Term
    )
  )

################################################################################

# Add Ontology definitions to second tab in spreadsheet.

# Find all unique ontology terms.
ontology_terms = bind_rows(
  dplyr::select(all, Basic_Data_Type_Ontology_Term) %>%
    rename(Term = Basic_Data_Type_Ontology_Term),
  dplyr::select(all, Field_Name_Ontology_Terms) %>%
    rename(Term = Field_Name_Ontology_Terms),
  dplyr::select(all, Value_Ontology_Terms) %>%
    rename(Term = Value_Ontology_Terms)) %>%
  separate_longer_delim(Term, delim = "||") %>%
  distinct(Term) %>%
  arrange(Term) %>%
  pull(Term)

# Read from ontology (version 26.04d of NCIT).
tmp_file_path = paste0(tempdir(), "/NCIT.csv.gz")
if (!file.exists(tmp_file_path)) {
  download.file("https://data.bioontology.org/ontologies/NCIT/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb&download_format=csv", tmp_file_path)
}
ncit = read_csv(tmp_file_path) %>%
  dplyr::select(code, `Preferred Label`, Synonyms, Definitions) %>%
  dplyr::rename(Code = code, Preferred_Label = `Preferred Label`) %>%
  mutate(Term = str_c(Preferred_Label, " (", Code, ")"))

# Did we get everything?
ncit_terms = pull(ncit, Term)
missing = setdiff(ontology_terms, ncit_terms)

if (length(missing) > 0) {
  print("These terms are not in the ontology.")
  print(ncit_terms)
  stop()
}

ncit = filter(ncit, Term %in% ontology_terms) %>%
  dplyr::select(-Term) %>%
  arrange(Preferred_Label)

write_xlsx(x = list(Mappings = all, Ontology_Terms = ncit), path = out_file_path)