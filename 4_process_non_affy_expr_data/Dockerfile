FROM bioconductor/bioconductor_docker:RELEASE_3_17

RUN R -e 'BiocManager::install(c("tidyverse", "GEOquery", "janitor"))'

WORKDIR /4_process_non_affy_expr_data
