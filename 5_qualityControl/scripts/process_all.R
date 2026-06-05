library(tidyverse)
library(tools)
library(doppelgangR)
library(ggplot2)
library(affy)
library(methods)
library(AnnotationDbi)
library(Biobase)
library(oligo)
library(GEOquery)

library(pd.hugene.1.0.st.v1)
library(pd.huex.1.0.st.v2)
library(u133aaofav2cdf)
library(hgu95av2cdf)
library(hgu133acdf)
library(hgu133a2cdf)
library(hgu133plus2cdf)
library(hgu133bcdf)

# This setting helps with the process very large files.
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10000)

options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

if (!dir.exists("/Data/expression_data2/")) {
  dir.create("/Data/expression_data2/")
}

if (!dir.exists("/Data/expression_data3/")) {
  dir.create("/Data/expression_data3/")
}

if (!dir.exists("/Data/expression_data4/")) {
  dir.create("/Data/expression_data4/")
}

if (!dir.exists("/Data/prelim_metadata2/")) {
  dir.create("/Data/prelim_metadata2/")
}

if (!dir.exists("/Data/analysis_ready_metadata")) {
  dir.create("/Data/analysis_ready_metadata")
}

if (!dir.exists("/Data/analysis_ready_expression_data")) {
  dir.create("/Data/analysis_ready_expression_data")
}

if (!dir.exists("/Data/doppelgangR_results")) {
  dir.create("/Data/doppelgangR_results")
}

if (!dir.exists("/Data/IQRray_results")) {
  dir.create("/Data/IQRray_results")
}

source("functions/compute_IQRray.R")
source("functions/run_IQRray.R")
source("functions/bind_IQR_file.R")

#source("scripts/clean_expression_data_colnames.R")
#source("scripts/match_expr_data_and_metadata.R")
#source("scripts/add_gene_identifiers.R")

# Separates datasets by array type
source("scripts/filter_chips.R")

source("scripts/doppelgangR.R")
###source("scripts/merge_doppel_results.R")

#source("scripts/IQRray_E_TABM_158.R")
source("scripts/IQRray_single_chips.R")
#source("scripts/IQRray_multiple_chips.R")
