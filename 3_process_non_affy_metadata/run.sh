#! /bin/bash

set -o errexit

#######################################################
# Build the Docker image
#######################################################

docker build -t inwosu/bc_data_curation_03 .

#######################################################
# Run docker command
#######################################################

dockerCommand="docker run -i -t --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/3_process_non_affy_metadata \
    -v $(pwd)/../Data:/Data \
    inwosu/bc_data_curation_03"

time $dockerCommand Rscript scripts/process_all.R

# $dockerCommand bash
