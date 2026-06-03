#! /bin/bash

set -o errexit

#######################################################
# Build the Docker image
#######################################################

docker build -t inwosu_bc_data_paper_map_metadata .

#######################################################
# Run docker command
#######################################################

# While testing, use this command:
dockerCommand="docker run -i -t --rm \
    -u $(id -u):$(id -g) \
    -v $(pwd):/7_map_metadata \
    -v $(pwd)/../Data:/Data \
    inwosu_bc_data_paper_map_metadata"

time $dockerCommand Rscript Check_All.R

# $dockerCommand bash
