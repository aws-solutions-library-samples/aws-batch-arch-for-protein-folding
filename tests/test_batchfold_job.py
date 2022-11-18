# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

import pytest
import boto3
from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.batchfold_job import BatchFoldJob
from time import sleep

@pytest.fixture()
def batch_environment():
    stack = BatchFoldEnvironment(boto_session=boto3.Session())
    return stack

def test_batchfold_job():

    job = BatchFoldJob()
    job.define_container_overrides(["echo hello"], cpu=8, memory=9, gpu=3)
    assert job.container_overrides["command"] == ["echo hello"]
    assert job.container_overrides["resourceRequirements"][0]["type"] == "VCPU"
    assert job.container_overrides["resourceRequirements"][0]["value"] == "8"
    assert job.container_overrides["resourceRequirements"][1]["type"] == "MEMORY"
    assert job.container_overrides["resourceRequirements"][1]["value"] == "9000"
    assert job.container_overrides["resourceRequirements"][2]["type"] == "GPU"
    assert job.container_overrides["resourceRequirements"][2]["value"] == "3"    
