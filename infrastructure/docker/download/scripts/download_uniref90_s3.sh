#!/bin/bash
#
# Original Copyright 2021 DeepMind Technologies Limited
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
#
# Downloads and unzips the UniRef90 database for AlphaFold.
#
# Usage: bash download_uniref90.sh /path/to/download/directory
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
ROOT_DIR="${DOWNLOAD_DIR}/uniref90"
SOURCE_URL="s3://aws-batch-architecture-for-alphafold-public-artifacts/uniref90/uniref90.fasta.gz"
BASENAME=$(basename "${SOURCE_URL}")

mkdir --parents "${ROOT_DIR}"
aws s3 cp --no-sign-request "${SOURCE_URL}" "${ROOT_DIR}"

pushd "${ROOT_DIR}"
gunzip "${ROOT_DIR}/${BASENAME}"
popd
