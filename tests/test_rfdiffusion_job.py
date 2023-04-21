# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.rfdiffusion_job import RFDiffusionJob
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

def test_rfdiffusion_job_init():
    bucket = os.getenv("TEST_BUCKET")
    job_name = "RFDiffusionJob" + datetime.now().strftime("%Y%m%d%s")
    params = {
        "inference.num_designs": "3",
        "contigmap.contigs": "[10-40/A163-181/10-40]"
    }
    new_job = RFDiffusionJob(
        boto_session=boto_session,
        job_name = job_name,
        memory=15,
        input_s3_uri = f"s3://{bucket}/rfdiffusion/5TPN/5TPN.pdb",
        output_s3_uri = f"s3://{bucket}/rfdiffusion/5TPN/outputs/",
        params=params
    )
    assert new_job.job_definition_name == "RFDiffusionJobDefinition"
    assert new_job.input_s3_uri == f"s3://{bucket}/rfdiffusion/5TPN/5TPN.pdb"
    assert new_job.output_s3_uri == f"s3://{bucket}/rfdiffusion/5TPN/outputs/"
    assert new_job.params["inference.num_designs"] == "3"
    assert new_job.params["contigmap.contigs"] == "[10-40/A163-181/10-40]"

def test_rfdiffusion_job_submission(batch_environment):
    bucket = os.getenv("TEST_BUCKET")
    job_queue_name = "G4dnJobQueue"
    job_name = "RFDiffusionJob" + datetime.now().strftime("%Y%m%d%s")
    params = {
        "inference.num_designs": "3",
        "contigmap.contigs": "[10-40/A163-181/10-40]"
    }
    new_job = RFDiffusionJob(
        boto_session=boto_session,
        job_name = job_name,
        memory=15,
        input_s3_uri = f"s3://{bucket}/rfdiffusion/5TPN/5TPN.pdb",
        output_s3_uri = f"s3://{bucket}/rfdiffusion/5TPN/outputs/",
        params = params
    )
    assert new_job.job_definition_name == "RFDiffusionJobDefinition"
    submission = batch_environment.submit_job(new_job, job_queue_name)
    assert job_name == submission.job_name    
    job_description = new_job.describe_job()        
    assert job_name == job_description[0].get("jobName", [])