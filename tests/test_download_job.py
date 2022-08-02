# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.download_job import DownloadJob
import pytest
import boto3
from datetime import datetime

@pytest.fixture()
def batch_environment():
    stack = BatchFoldEnvironment(boto_session = boto3.Session())
    assert "BatchEnvironment" in stack.stack_name
    return(stack)

def test_download_job_init():
    new_job = DownloadJob(
        job_definition_name = "Test",
        script = "test.sh",
        data_dir = "/test"
    )

    assert new_job.job_definition_name == "Test"
    assert new_job.script == "test.sh"
    assert new_job.data_dir == "/test"

def test_download_job_submission(batch_environment):
    job_name = "DownloadJob" + datetime.now().strftime("%Y%m%d%s")
    job_queue_name = "GravitonSpotJobQueue"

    new_job = DownloadJob(
        job_name = job_name, 
        script = "./scripts/download_test.sh",
        cpu = 4
    )

    submission = batch_environment.submit_job(new_job, job_queue_name)
    assert job_name == submission.job_name
    
    job_description = new_job.describe_job()        
    assert job_name == job_description[0].get("jobName", [])

    job_dict = batch_environment.list_jobs()
    job_list = job_dict[job_queue_name]
    assert len(job_list) > 0

    job_info = [job for job in job_list if job.get("jobName", []) == job_name]
    assert job_info[0].get("jobDefinition") == batch_environment.job_definitions["MSAJobDefinition"]
    assert job_info[0].get("jobName") == job_name