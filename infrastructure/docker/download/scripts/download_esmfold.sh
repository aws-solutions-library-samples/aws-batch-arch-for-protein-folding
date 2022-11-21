#!/bin/bash
#
# Original Copyright 2022 AlQuraishi Laboratory
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
#
# Downloads and unzips the RFDesign parameters.
#
# Usage: bash download_esmfold.sh /path/to/download/directory
set -e

if [[ $# -eq 0 ]]; then
    echo "Error: download directory must be provided as an input argument."
    exit 1
fi

if ! command -v aws &> /dev/null ; then
    echo "Error: awscli could not be found. Please install awscli."
    exit 1
fi

DOWNLOAD_DIR="${1}/esmfold_params/hub/checkpoints"
mkdir -p "${DOWNLOAD_DIR}"
aria2c "https://dl.fbaipublicfiles.com/fair-esm/models/esmfold_3B_v1.pt" --dir=${DOWNLOAD_DIR}
aria2c "https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t36_3B_UR50D.pt" --dir=${DOWNLOAD_DIR}
aria2c "https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t36_3B_UR50D-contact-regression.pt" --dir=${DOWNLOAD_DIR}
