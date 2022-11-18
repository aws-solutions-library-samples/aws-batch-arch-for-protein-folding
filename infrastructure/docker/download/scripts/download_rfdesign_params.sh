#!/bin/bash
#
# Original Copyright 2022 AlQuraishi Laboratory
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
#
# Downloads and unzips the RFDesign parameters.
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

DOWNLOAD_DIR="${1}/rfdesign_params"
mkdir -p "${DOWNLOAD_DIR}"
aria2c "http://files.ipd.uw.edu/pub/rfdesign/weights/BFF_last.pt" --dir="${DOWNLOAD_DIR}/hallucination/rf_Nov05"
aria2c "https://raw.githubusercontent.com/RosettaCommons/RFDesign/main/hallucination/weights/rf_Nov05/params_Nov05.json" --dir="${DOWNLOAD_DIR}/hallucination/rf_Nov05"
aria2c "http://files.ipd.uw.edu/pub/rfdesign/weights/BFF_mix_epoch25.pt" --dir="${DOWNLOAD_DIR}/inpainting"