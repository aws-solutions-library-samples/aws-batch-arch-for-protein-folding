#!/bin/bash
#
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
#
# Downloads and prepares MMSeqs2-formatted sequence databases.
#
# This script downloads UniRef90 and SwissProt databases in MMSeqs2's native
# indexed format using the `mmseqs databases` command.
#
# Note: BFD and UniRef30 databases require additional conversion steps.
# For these databases, use the ColabFold setup script available at:
# https://colabfold.mmseqs.com/
#
# Usage: bash download_mmseqs2_dbs.sh /path/to/download/directory
set -e

if [[ $# -eq 0 ]]; then
    echo "Error: download directory must be provided as an input argument."
    exit 1
fi

DOWNLOAD_DIR="$1"
ROOT_DIR="${DOWNLOAD_DIR}/mmseqs2_dbs"

mkdir --parents "${ROOT_DIR}"

# Install mmseqs if not already present
if ! command -v mmseqs &> /dev/null; then
    echo "Installing MMSeqs2..."
    wget -q https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz -O /tmp/mmseqs.tar.gz
    tar xzf /tmp/mmseqs.tar.gz -C /tmp
    mv /tmp/mmseqs/bin/mmseqs /usr/local/bin/
    chmod +x /usr/local/bin/mmseqs
    rm -rf /tmp/mmseqs /tmp/mmseqs.tar.gz
fi

# Download UniRef90 in MMSeqs2 format
echo "Downloading UniRef90 in MMSeqs2 format..."
UNIREF90_DIR="${ROOT_DIR}/uniref90"
mkdir --parents "${UNIREF90_DIR}"
cd "${UNIREF90_DIR}"
if [ ! -f "uniref90.dbtype" ]; then
    mmseqs databases UniRef90 uniref90 tmp --threads 4
    rm -rf tmp
    echo "UniRef90 download complete."
else
    echo "UniRef90 already exists, skipping."
fi

# Download SwissProt in MMSeqs2 format
echo "Downloading SwissProt in MMSeqs2 format..."
SWISSPROT_DIR="${ROOT_DIR}/swissprot"
mkdir --parents "${SWISSPROT_DIR}"
cd "${SWISSPROT_DIR}"
if [ ! -f "swissprot.dbtype" ]; then
    mmseqs databases "UniProtKB/Swiss-Prot" swissprot tmp --threads 4
    rm -rf tmp
    echo "SwissProt download complete."
else
    echo "SwissProt already exists, skipping."
fi

echo "MMSeqs2 database downloads complete."
echo "Downloaded databases are located at: ${ROOT_DIR}"
echo ""
echo "Note: For BFD and UniRef30 databases in MMSeqs2 format, use the"
echo "ColabFold setup script: https://colabfold.mmseqs.com/"
