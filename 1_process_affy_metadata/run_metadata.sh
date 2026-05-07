#! /bin/bash

set -o errexit

#######################################################
# Build the Docker image
#######################################################

docker build -t inwosu/bc_data_curation_01 .

#######################################################
# Run docker command
#######################################################
 
dockerCommand="docker run -i -t --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/1_process_metadata \
    -v $(pwd)/../Data:/Data \
    inwosu/bc_data_curation_01"

time $dockerCommand Rscript scripts/parse_metadata.R 

# $dockerCommand bash