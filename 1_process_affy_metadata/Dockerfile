FROM bioconductor/bioconductor_docker:RELEASE_3_18

#install required packages
RUN R -e 'BiocManager::install(c("GEOquery", "tidyverse", "stringi", "janitor", "rlist"), force = TRUE)'

ADD *.tsv /

WORKDIR /1_process_metadata
