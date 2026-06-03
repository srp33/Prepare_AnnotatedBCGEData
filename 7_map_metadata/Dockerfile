FROM bioconductor/bioconductor_docker:RELEASE_3_17

RUN R -e 'BiocManager::install(c("tidyverse", "writexl"), force = TRUE)'

WORKDIR /7_map_metadata
