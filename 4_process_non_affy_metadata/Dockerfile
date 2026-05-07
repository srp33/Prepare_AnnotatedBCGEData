FROM bioconductor/bioconductor_docker:RELEASE_3_17

RUN R -e 'BiocManager::install(c("GEOquery", "tidyverse", "janitor"))'

WORKDIR /4_process_non_affy_metadata
