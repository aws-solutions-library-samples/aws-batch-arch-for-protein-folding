# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.proteinmpnn_job import ProteinMPNNJob
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

def test_proteinmpnn_job_init():

    bucket = os.getenv("TEST_BUCKET")
    job_name = "ProteinMPNNJob" + datetime.now().strftime("%Y%m%d%s")

    new_job = ProteinMPNNJob(
        boto_session=boto_session,
        job_name = job_name,
        memory=15,
        pdb_s3_uri = f"s3://{bucket}/mpnn/input/3HTN.pdb",
        output_s3_uri = f"s3://{bucket}/mpnn/output",
        pdb_path_chains = "A B",
        num_seq_per_target = 2,
        sampling_temp = 0.01,
        seed = 37,
        batch_size = 1
    )
    assert new_job.job_definition_name == "ProteinMPNNJobDefinition"
    assert new_job.pdb_s3_uri == f"s3://{bucket}/mpnn/input/3HTN.pdb"
    assert new_job.output_s3_uri == f"s3://{bucket}/mpnn/output"
    assert new_job.pdb_path_chains == "A B"
    assert new_job.num_seq_per_target == 2
    assert new_job.sampling_temp == 0.01
    assert new_job.seed == 37
    assert new_job.batch_size == 1

def test_proteinmpnn_job_submission(batch_environment):

    job_name = "ProteinMPNNJob" + datetime.now().strftime("%Y%m%d%s")
    job_queue_name = "CPUOnDemandJobQueue"
    bucket = os.getenv("TEST_BUCKET")
    new_job = ProteinMPNNJob(
        boto_session=boto_session,
        job_name = job_name,
        memory=15,
        pdb_s3_uri = f"s3://{bucket}/mpnn/input/3HTN.pdb",
        output_s3_uri = f"s3://{bucket}/mpnn/output/{job_name}",
        pdb_path_chains = "A B",
        num_seq_per_target = 25,
        sampling_temp = 0.01,
        seed = 37,
        batch_size = 1
    )
    submission = batch_environment.submit_job(new_job, job_queue_name)
    assert job_name == submission.job_name    
    job_description = new_job.describe_job()        
    assert job_name == job_description[0].get("jobName", [])
    job_info = []
    while job_info == []:
        sleep(5)
        job_dict = batch_environment.list_jobs()
        job_list = job_dict[job_queue_name]  
        job_info = [job for job in job_list if job.get("jobName", []) == job_name]
    assert job_info[0].get("jobDefinition") == batch_environment.job_definitions["ProteinMPNNJobDefinition"]
    assert job_info[0].get("jobName") == job_name