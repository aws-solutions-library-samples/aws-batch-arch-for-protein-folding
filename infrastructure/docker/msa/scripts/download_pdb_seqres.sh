#!/bin/bash
#
# Original Copyright 2021 DeepMind Technologies Limited
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
#
# Downloads and unzips the PDB SeqRes database for AlphaFold.
#
# Usage: bash download_pdb_seqres.sh /path/to/download/directory
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
ROOT_DIR="${DOWNLOAD_DIR}/pdb_seqres"
LATEST_PDB_SNAPSHOT=$(aws s3 ls --no-sign-request s3://pdbsnapshots/ | tail -n 3 | head -n 1 | awk '{print $2}')

SOURCE_URL="s3://pdbsnapshots/${LATEST_PDB_SNAPSHOT}pub/pdb/derived_data/pdb_seqres.txt"
BASENAME=$(basename "${SOURCE_URL}")
mkdir --parents "${ROOT_DIR}"

aws s3 cp --no-sign-request "${SOURCE_URL}" "${ROOT_DIR}"
