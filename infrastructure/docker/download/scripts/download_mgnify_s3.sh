#!/bin/bash
#
# Original Copyright 2021 DeepMind Technologies Limited
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
#
# Downloads and unzips the MGnify database for AlphaFold.
#
# Usage: bash download_mgnify.sh /path/to/download/directory
set -e

if [[ $# -eq 0 ]]; then
    echo "Error: download directory must be provided as an input argument."
    exit 1
fi

if ! command -v aws &> /dev/null ; then
    echo "Error: awscli could not be found. Please install awscli."
    exit 1
fi

DOWNLOAD_DIR="$1"
ROOT_DIR="${DOWNLOAD_DIR}/mgnify"
# Mirror of:
# ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2018_12/mgy_clusters.fa.gz
SOURCE_URL="s3://aws-batch-architecture-for-alphafold-public-artifacts/mgnify/mgy_clusters_2018_12.fa.gz"
BASENAME=$(basename "${SOURCE_URL}")

mkdir --parents "${ROOT_DIR}"
aws s3 cp --no-sign-request "${SOURCE_URL}" "${ROOT_DIR}"
pushd "${ROOT_DIR}"
gunzip "${ROOT_DIR}/${BASENAME}"
popd
