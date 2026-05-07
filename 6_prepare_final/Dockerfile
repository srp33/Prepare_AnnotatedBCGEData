FROM bioconductor/bioconductor_docker:RELEASE_3_15

#install required packages
RUN R -e 'BiocManager::install(c("tidyverse", "janitor", "biomaRt"), force = TRUE)'

WORKDIR /6_prepare_final