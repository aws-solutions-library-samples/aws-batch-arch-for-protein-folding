# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.nextflow_job import NextFlowJob
import pytest
import boto3
from datetime import datetime
import os
from time import sleep

boto_session = boto3.Session(region_name=os.getenv("AWS_REGION"))

@pytest.fixture()
def batch_environment():
    stack = BatchFoldEnvironment(boto_session = boto_session)
    return(stack)

def test_nextflow_job_init():
    bucket = os.getenv("TEST_BUCKET")
    job_name = "NextFlowJob" + datetime.now().strftime("%Y%m%d%s")
    params = {
        "parameterA": "A",
        "parameterB": "B"
    }
    new_job = NextFlowJob(
        boto_session=boto_session,
        job_name = job_name,
        assets_s3_uri = os.path.join("s3://", bucket, "assets"),
        nextflow_script = 'test.nf',
        params = params
    )
    assert new_job.job_definition_name == "NextflowJobDefinition"
    assert new_job.assets_s3_uri == os.path.join("s3://", bucket, "assets")
    assert new_job.nextflow_script == 'test.nf'
    assert new_job.params["parameterA"] == "A"
    assert new_job.params["parameterB"] == "B"