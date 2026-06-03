#! /bin/bash

set -o errexit

#######################################################
# Build the Docker image
#######################################################

docker build -t inwosu/bc_data_curation_05 .

#######################################################
# Run docker command
#######################################################

# While you are testing, use this command:
#dockerCommand="docker run -i -t --rm \
dockerCommand="docker run -i --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/5_process_non_affy_expr_data \
    -v $(pwd)/../Data:/Data \
    inwosu/bc_data_curation_05"

time $dockerCommand Rscript scripts/process_all.R

# $dockerCommand bash
