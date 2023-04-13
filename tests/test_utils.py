# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from batchfold.utils import utils
import pytest
import boto3
from datetime import datetime
import os
from time import sleep

def test_download_rcsb_pdb_file():
    filename = utils.download_rcsb_pdb_file("5TPN", "tests/data", file_format="pdb")
    assert filename == "tests/data/5TPN.pdb"

def test_extract_pdb_chain():
    input_file = "tests/data/5TPN.pdb"
    model = 0
    chain = "H"
    filename = utils.extract_chain(input_file, model, chain)
    assert filename == "tests/data/5TPN_H.pdb"
