#! /bin/bash

set -o errexit

#######################################################
# Build the Docker image
#######################################################

cp ../*.tsv .
docker build -t inwosu/bc_data_curation_02 .  

#######################################################
# Run detailed functional tests on small file
#######################################################

#dockerCommand="docker run -i -t --rm \
dockerCommand="docker run -i --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/2_process_affy_expression_data \
    -v $(pwd)/../Data:/Data \
    inwosu/bc_data_curation_02"

time $dockerCommand Rscript scripts/parse_nomalize_scripts.R
 
# $dockerCommand bash
