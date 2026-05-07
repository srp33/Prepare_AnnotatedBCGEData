#! /bin/bash

set -o errexit

#######################################################
# Build the Docker image
#######################################################

docker build -t inwosu/bc_data_curation_04 .

#######################################################
# Run docker command
#######################################################

dockerCommand="docker run -i -t --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/4_process_non_affy_metadata \
    -v $(pwd)/../Data:/Data \
    inwosu/bc_data_curation_04"

time $dockerCommand Rscript scripts/source_all_non_affy_meta.R

# $dockerCommand bash
