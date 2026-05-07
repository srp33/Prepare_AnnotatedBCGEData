FROM bioconductor/bioconductor_docker:RELEASE_3_17

#install required packages
RUN R -e 'BiocManager::install(c("tidyverse", "doppelgangR", "affy", "AnnotationDbi", "Biobase", "oligo", "GEOquery", "ggplot2"), force = TRUE)'
RUN R -e 'BiocManager::install(c("pd.hugene.1.0.st.v1", "pd.huex.1.0.st.v2", "u133aaofav2cdf", "hgu95av2cdf", "hgu133acdf", "hgu133a2cdf", "hgu133plus2cdf", "hgu133bcdf"), force = TRUE)'

WORKDIR /3_qualityControl_affy_expression_data


# pd.hugene.1.0.st.v1   Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]
# pd.huex.1.0.st.v2     Affymetrix Human Exon 1.0 ST Array [transcript (gene) version]
# u133aaofav2cdf        Affymetrix GeneChip HT-HG_U133A Early Access Array
# hgu95av2cdf           Affymetrix Human Genome U95 Version 2 Array
# hgu133acdf            Affymetrix Human Genome U133A Array
# hgu133a2cdf           Affymetrix Human Genome U133A 2.0 Array
# hgu133plus2cdf        Affymetrix Human Genome U133 Plus 2.0 Array
# hgu133bcdf            Affymetrix Human Genome U133B Array