#!/bin/bash
#
# Original Copyright 2021 DeepMind Technologies Limited
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
#
# Downloads, unzips and flattens the PDB database for AlphaFold.
#
# Usage: bash download_pdb_mmcif.sh /path/to/download/directory
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
ROOT_DIR="${DOWNLOAD_DIR}/pdb_mmcif"
RAW_DIR="${ROOT_DIR}/raw"
MMCIF_DIR="${ROOT_DIR}/mmcif_files"

mkdir --parents "${RAW_DIR}"

LATEST_PDB_SNAPSHOT=$(aws s3 ls --no-sign-request s3://pdbsnapshots/ | tail -n 3 | head -n 1 | awk '{print $2}')

aws s3 cp --recursive --no-sign-request s3://pdbsnapshots/"${LATEST_PDB_SNAPSHOT}"pub/pdb/data/structures/divided/mmCIF "${RAW_DIR}"

echo "Unzipping all mmCIF files..."
find "${RAW_DIR}/" -type f -iname "*.gz" -exec gunzip {} +

echo "Flattening all mmCIF files..."
mkdir --parents "${MMCIF_DIR}"
find "${RAW_DIR}" -type d -empty -delete  # Delete empty directories.
for subdir in "${RAW_DIR}"/*; do
  mv "${subdir}/"*.cif "${MMCIF_DIR}"
done

# Delete empty download directory structure.
find "${RAW_DIR}" -type d -empty -delete

aws s3 cp --no-sign-request s3://pdbsnapshots/"${LATEST_PDB_SNAPSHOT}"pub/pdb/data/status/obsolete.dat "${ROOT_DIR}"