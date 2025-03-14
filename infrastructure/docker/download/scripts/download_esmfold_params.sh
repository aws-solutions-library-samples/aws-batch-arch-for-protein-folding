#!/bin/bash
#
# Original Copyright 2022 AlQuraishi Laboratory
# Modifications Copyright 2023 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
#
# Downloads and unzips the ESMFold parameters.
#
# Usage: bash download_esmfold.sh /path/to/download/directory
set -e

if [[ $# -eq 0 ]]; then
    echo "Error: download directory must be provided as an input argument."
    exit 1
fi

DOWNLOAD_DIR="$1"
ROOT_DIR="${DOWNLOAD_DIR}/esmfold_params/hub"
SOURCE_URL="facebook/esmfold_v1"

git lfs install
echo "Downloading ${SOURCE_URI} from Hugging Face Hub"
mkdir tmp
git clone https://huggingface.co/$SOURCE_URI tmp --depth=1
rm -rf tmp/.git
mv -n tmp/* $ROOT_DIR
rm -rf tmp
