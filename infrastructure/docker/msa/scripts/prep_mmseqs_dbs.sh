#!/bin/bash
#
# Original Copyright 2022 AlQuraishi Laboratory
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
#
#Processes the mmseqs databases

set -e

DOWNLOAD_DIR="$1"
ROOT_DIR="${DOWNLOAD_DIR}/mmseqs_dbs"
mkdir -p $ROOT_DIR

for f in $(ls ${DOWNLOAD_DIR}/*.tar.gz)
do
  tar --extract --verbose --file="${f}" \
      --directory=$ROOT_DIR
  rm "${f}"
  BASENAME="$(basename ${f%%.*})"
  DB_NAME="${BASENAME}_db"
  OLD_PWD=$(pwd)
  cd $ROOT_DIR 
  mmseqs tsv2exprofiledb "${BASENAME}" "${DB_NAME}"
  mmseqs createindex "${DB_NAME}" "${DOWNLOAD_DIR}/tmp/"
  cd "${OLD_PWD}"
done