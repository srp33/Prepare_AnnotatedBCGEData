FROM bioconductor/bioconductor_docker:RELEASE_3_17

#install required packages
RUN R -e 'BiocManager::install(c("GEOquery", "SCAN.UPC", "tidyverse", "janitor", "stringi", "biomaRt", "rlist", "doParallel", "doppelgangR", "tools"), force = TRUE)'

#add brain array script and install packages
ADD brain_Array.R / 
RUN Rscript /brain_Array.R

WORKDIR /2_process_affy_expression_data
