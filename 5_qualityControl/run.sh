#! /bin/bash

set -o errexit

#######################################################
# Build the Docker image
#######################################################

docker build -t inwosu/bc_data_curation_05 .

#######################################################
# Run detailed functional tests on small file
#######################################################

#dockerCommand="docker run -i -t --rm \
dockerCommand="docker run -i --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/5_qualityControl \
    -v $(pwd)/../Data:/Data \
    inwosu/bc_data_curation_05"

time $dockerCommand Rscript scripts/process_all.R

# $dockerCommand bash    
