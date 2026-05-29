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
# Add ontology terms for the primitive type of each (categorical, numeric, maybe not integer or identifier).

# Manual check.
#   Pivot longer and then evaluate whether the same fields and values seem to
#   have the same ontology terms mapped to them.
# Add the unique values for numeric and identifier. Pivot longer.