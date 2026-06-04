# Identify and remove duplicate rows.
source("Remove_Duplicates.R")

# Find fields that are annotated in both categorical and numeric. In some cases,
#   this is fine because there are some categorical values and some numeric.
source("Find_Problems.R")

# Identify which columns are in the metadata but have not been mapped.
source("Check_Fields.R")

# Identify which categorical values are in the metadata but have not been mapped.
source("Check_Categorical_Values.R")

# Make sure the codes have been specified properly in the spreadsheet.
source("Check_Annotations_Syntax.R")

# Make sure the ontology preferred terms and codes align properly.
source("Check_Ontology_Alignment.R")

# Peform miscellaneous cleaning tasks.
source("Clean.R")

# Manual check.
#   Sort and review the final spreadsheet various ways.
#   Evaluate whether the same fields and values
#   have the same ontology terms mapped to them, etc.