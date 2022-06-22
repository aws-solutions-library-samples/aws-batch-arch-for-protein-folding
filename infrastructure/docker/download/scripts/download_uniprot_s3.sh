#!/bin/bash
#
# Original Copyright 2021 DeepMind Technologies Limited
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
#
# Downloads, unzips and merges the SwissProt and TrEMBL databases for
# AlphaFold-Multimer.
#
# Usage: bash download_uniprot.sh /path/to/download/directory
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
ROOT_DIR="${DOWNLOAD_DIR}/uniprot"

TREMBL_SOURCE_URL="s3://aws-batch-architecture-for-alphafold-public-artifacts/uniprot/uniprot_trembl.fasta.gz"
TREMBL_BASENAME=$(basename "${TREMBL_SOURCE_URL}")
TREMBL_UNZIPPED_BASENAME="${TREMBL_BASENAME%.gz}"

SPROT_SOURCE_URL="s3://aws-batch-architecture-for-alphafold-public-artifacts/uniprot/uniprot_sprot.fasta.gz"
SPROT_BASENAME=$(basename "${SPROT_SOURCE_URL}")
SPROT_UNZIPPED_BASENAME="${SPROT_BASENAME%.gz}"

mkdir --parents "${ROOT_DIR}"
aws s3 cp --no-sign-request "${TREMBL_SOURCE_URL}" "${ROOT_DIR}"
aws s3 cp --no-sign-request "${SPROT_SOURCE_URL}" "${ROOT_DIR}"

pushd "${ROOT_DIR}"
gunzip "${ROOT_DIR}/${TREMBL_BASENAME}"
gunzip "${ROOT_DIR}/${SPROT_BASENAME}"

# Concatenate TrEMBL and SwissProt, rename to uniprot and clean up.
cat "${ROOT_DIR}/${SPROT_UNZIPPED_BASENAME}" >> "${ROOT_DIR}/${TREMBL_UNZIPPED_BASENAME}"
mv "${ROOT_DIR}/${TREMBL_UNZIPPED_BASENAME}" "${ROOT_DIR}/uniprot.fasta"
rm "${ROOT_DIR}/${SPROT_UNZIPPED_BASENAME}"
popd
