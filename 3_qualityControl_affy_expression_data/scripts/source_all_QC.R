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

library("pd.hugene.1.0.st.v1")
library("pd.huex.1.0.st.v2")
library("u133aaofav2cdf")
library("hgu95av2cdf")
library("hgu133acdf")
library("hgu133a2cdf")
library("hgu133plus2cdf")
library("hgu133bcdf")

options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

doppel_dir <- "/Data/doppelgang_results/"
    if (!dir.exists(doppel_dir)) {
      dir.create(doppel_dir)
    }

IQRray_file_path <- "/Data/IQRray_results/"
if (!dir.exists(IQRray_file_path)) {
  dir.create(IQRray_file_path)
}

# source required functions
source("functions/compute_IQRray.R")
source("functions/run_IQRray.R")
source("functions/bind_IQR_file.R")

#source script that seperates datasets by array type
source("scripts/filter_chips.R")

# various QC scripts
source("scripts/doppelgang.R")
source("scripts/merge_doppel_results.R")
source("scripts/IQRray_E_TABM_158.R")
source("scripts/IQRray_single_chips.R")
source("scripts/IQRray_multiple_chips.R")

# delete temporary download directory
unlink("GSE1456", recursive = TRUE, force = TRUE)
unlink("GSE3494", recursive = TRUE, force = TRUE)
unlink("GSE4922", recursive = TRUE, force = TRUE)
unlink("GSE6532", recursive = TRUE, force = TRUE)