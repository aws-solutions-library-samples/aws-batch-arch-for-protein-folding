# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.diffdock_job import DiffDockJob
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

def test_diffdock_job_init():
    bucket = os.getenv("TEST_BUCKET")
    new_job = DiffDockJob(
        boto_session = boto_session,
        protein_s3_uri = f"s3://{bucket}/6W70/6w70.pdb",
        ligand_s3_uri = f"s3://{bucket}/6W70/6w70_ligand.sdf",
        output_s3_uri = f"s3://{bucket}/6W70/outputs/",
        complex_name = "6W70",
    )
    assert new_job.job_definition_name == "DiffDockJobDefinition"
    assert new_job.protein_s3_uri == f"s3://{bucket}/6W70/6w70.pdb"
    assert new_job.ligand_s3_uri == f"s3://{bucket}/6W70/6w70_ligand.sdf"
    assert new_job.output_s3_uri == f"s3://{bucket}/6W70/outputs/"
    assert new_job.complex_name == "6W70"
    assert new_job.save_visualisation
    assert new_job.samples_per_complex == 10

def test_diffdock_job_submission(batch_environment):

    job_name = "DiffDockJob" + datetime.now().strftime("%Y%m%d%s")
    job_queue_name = "G4dnJobQueue"
    bucket = os.getenv("TEST_BUCKET")
    new_job = DiffDockJob(
        job_name = job_name,
        protein_s3_uri = f"s3://{bucket}/6W70/6w70.pdb",
        ligand_s3_uri = f"s3://{bucket}/6W70/6w70_ligand.sdf",
        output_s3_uri = f"s3://{bucket}/6W70/outputs/",
        complex_name = "6W70",
    )
    submission = batch_environment.submit_job(new_job, job_queue_name)
    assert job_name == submission.job_name    
    job_description = new_job.describe_job()        
    # assert job_name == job_description[0].get("jobName", [])
    job_info = []
    while job_info == []:
        sleep(5)
        job_dict = batch_environment.list_jobs()
        job_list = job_dict[job_queue_name]  
        job_info = [job for job in job_list if job.get("jobName", []) == job_name]
    assert job_info[0].get("jobDefinition") == batch_environment.job_definitions["DiffDockJobDefinition"]
    assert job_info[0].get("jobName") == job_name