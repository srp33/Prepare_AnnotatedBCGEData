FROM bioconductor/bioconductor_docker:RELEASE_3_23

# Install only packages that are not already present. Fail the build if any
# package cannot be loaded, so a failed install is not cached as success.
# Written with printf so older Docker versions, which do not support heredocs, can parse this file.
RUN printf '%s\n' \
  'install_pkgs <- function(pkgs) {' \
  '  missing_pkgs <- setdiff(pkgs, rownames(installed.packages()))' \
  '  if (length(missing_pkgs)) {' \
  '    BiocManager::install(missing_pkgs, update = FALSE, ask = FALSE)' \
  '  }' \
  '  failed <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]' \
  '  if (length(failed)) {' \
  '    message("Failed to install: ", paste(failed, collapse = ", "))' \
  '    quit(status = 1)' \
  '  }' \
  '}' \
  > /usr/local/bin/install_pkgs.R

RUN Rscript -e 'source("/usr/local/bin/install_pkgs.R"); install_pkgs(c("tidyverse", "data.table", "AnnotationDbi", "Biobase", "doppelgangR", "affy", "oligo", "GEOquery", "doParallel"))'
RUN Rscript -e 'source("/usr/local/bin/install_pkgs.R"); install_pkgs(c("pd.hugene.1.0.st.v1", "pd.huex.1.0.st.v2", "u133aaofav2cdf", "hgu95av2cdf", "hgu133acdf", "hgu133a2cdf", "hgu133plus2cdf", "hgu133bcdf"))'
RUN Rscript -e 'source("/usr/local/bin/install_pkgs.R"); install_pkgs("arrayQualityMetrics")'
RUN Rscript -e 'source("/usr/local/bin/install_pkgs.R"); install_pkgs("variancePartition")'

WORKDIR /5_qualityControl

# pd.hugene.1.0.st.v1   Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]
# pd.huex.1.0.st.v2     Affymetrix Human Exon 1.0 ST Array [transcript (gene) version]
# u133aaofav2cdf        Affymetrix GeneChip HT-HG_U133A Early Access Array
# hgu95av2cdf           Affymetrix Human Genome U95 Version 2 Array
# hgu133acdf            Affymetrix Human Genome U133A Array
# hgu133a2cdf           Affymetrix Human Genome U133A 2.0 Array
# hgu133plus2cdf        Affymetrix Human Genome U133 Plus 2.0 Array
# hgu133bcdf            Affymetrix Human Genome U133B Array
