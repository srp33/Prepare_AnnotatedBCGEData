#! /bin/bash

set -o errexit

#######################################################
# Build the Docker image
#######################################################

docker build -t inwosu/bc_data_curation_03 .

#######################################################
# Run detailed functional tests on small file
#######################################################

dockerCommand="docker run -i -t --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/3_qualityControl_affy_expression_data \
    -v $(pwd)/../Data:/Data \
    inwosu/bc_data_curation_03"

time $dockerCommand Rscript scripts/process_all.R

# $dockerCommand bash    
