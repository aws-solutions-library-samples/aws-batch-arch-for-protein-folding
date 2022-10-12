#!/bin/bash
#
# Original Copyright 2022 AlQuraishi Laboratory
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
#
# Downloads and unzips the OpenFold parameters.
#
# Usage: bash download_openfold_params.sh /path/to/download/directory
set -e

if [[ $# -eq 0 ]]; then
    echo "Error: download directory must be provided as an input argument."
    exit 1
fi

if ! command -v aws &> /dev/null ; then
    echo "Error: awscli could not be found. Please install awscli."
    exit 1
fi

DOWNLOAD_DIR="${1}/omegafold_params"
mkdir -p "${DOWNLOAD_DIR}"
aws s3 cp --no-sign-request s3://helixon/release1.pt "${DOWNLOAD_DIR}"/model.pt