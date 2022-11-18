# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.rfdesign_job import RFDesignHallucinateJob, RFDesignInpaintJob
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

def test_rfdesign_hallucinate_job_init():
    bucket = os.getenv("TEST_BUCKET")
    job_name = "RFDesignHallucinateJob" + datetime.now().strftime("%Y%m%d%s")
    mask = '25-35,B63-82,15-25,B119-140,0-15' 
    params = {
        "mask": mask,
        "steps": "g3",
        "w_rog": 1,
        "rog_thresh": 16,
        "w_rep": 2,
        "rep_pdb": "input/pdl1.pdb",
        "rep_sigma": 4,
        "save_pdb": False,
        "track_step": 10
    }
    new_job = RFDesignHallucinateJob(
        boto_session=boto_session,
        job_name = job_name,
        memory=15,
        target_id = "4ZQK",
        input_s3_uri = f"s3://{bucket}/pd1-demo/",
        output_s3_uri = f"s3://{bucket}/pd1-demo/outputs/",
        pdb = "input/pd1.pdb",
        params = params
    )
    assert new_job.job_definition_name == "RFDesignJobDefinition"
    assert new_job.target_id == "4ZQK"
    assert new_job.input_s3_uri == f"s3://{bucket}/pd1-demo/"
    assert new_job.output_s3_uri == f"s3://{bucket}/pd1-demo/outputs/"
    assert new_job.pdb == "input/pd1.pdb"
    assert new_job.params['mask'] == mask
    assert new_job.params['w_rep'] == 2
    assert new_job.params['save_pdb'] == False

def test_rfdesign_hallucinate_job_submission(batch_environment):

    job_name = "RFDesignHallucinateJob" + datetime.now().strftime("%Y%m%d%s")
    job_queue_name = "G4dnJobQueue"
    bucket = os.getenv("TEST_BUCKET")
    mask = '25-35,B63-82,15-25,B119-140,0-15' 
    params = {
        "mask": mask,
        "steps": "g3",
        "w_rog": 1,
        "rog_thresh": 16,
        "w_rep": 2,
        "rep_pdb": "input/pdl1.pdb",
        "rep_sigma": 4,
        "save_pdb": False,
        "track_step": 10
    }
    new_job = RFDesignHallucinateJob(
        boto_session=boto_session,
        job_name = job_name,
        memory=15,
        target_id = "4ZQK",
        input_s3_uri = f"s3://{bucket}/pd1-demo/",
        output_s3_uri = f"s3://{bucket}/pd1-demo/outputs/",
        pdb = "input/pd1.pdb",
        params = params
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
    assert job_info[0].get("jobDefinition") == batch_environment.job_definitions["RFDesignJobDefinition"]
    assert job_info[0].get("jobName") == job_name

def test_rfdesign_inpainting_job_init():
    bucket = os.getenv("TEST_BUCKET")
    job_name = "RFDesignInpaintingJob" + datetime.now().strftime("%Y%m%d%s")
    params = {
        "contigs":"25-35,B63-82,15-25,B119-140,0-15",
        "len": "80-115",
        "num_designs": 4,
        "dump_all": True,
    }
    new_job = RFDesignInpaintJob(
        boto_session=boto_session,
        job_name = job_name,
        memory=15,
        target_id = "4ZQK",
        input_s3_uri = f"s3://{bucket}/pd1-demo/pd1.pdb",
        output_s3_uri = f"s3://{bucket}/pd1-demo/outputs/",
        pdb = "input/pd1.pdb",
        params = params
    )
    assert new_job.job_definition_name == "RFDesignJobDefinition"
    assert new_job.target_id == "4ZQK"
    assert new_job.input_s3_uri == f"s3://{bucket}/pd1-demo/pd1.pdb"
    assert new_job.output_s3_uri == f"s3://{bucket}/pd1-demo/outputs/"
    assert new_job.pdb == "input/pd1.pdb"
    assert new_job.params['contigs'] == "25-35,B63-82,15-25,B119-140,0-15"
    assert new_job.params['len'] == "80-115"
    assert new_job.params['num_designs'] == 4
    assert new_job.params['dump_all'] == True

def test_rfdesign_inpainting_job_submission(batch_environment):

    job_name = "RFDesignInpaintingJob" + datetime.now().strftime("%Y%m%d%s")
    job_queue_name = "G4dnJobQueue"
    bucket = os.getenv("TEST_BUCKET")
    params = {
        "contigs":"25-35,B63-82,15-25,B119-140,0-15",
        "len": "80-115",
        "num_designs": 4,
        "dump_all": True,
    }
    new_job = RFDesignInpaintJob(
        boto_session=boto_session,
        job_name = job_name,
        memory=15,
        target_id = "4ZQK",
        input_s3_uri = f"s3://{bucket}/pd1-demo/pd1.pdb",
        output_s3_uri = f"s3://{bucket}/pd1-demo/outputs/",
        pdb = "input/pd1.pdb",
        params = params
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
    assert job_info[0].get("jobDefinition") == batch_environment.job_definitions["RFDesignJobDefinition"]
    assert job_info[0].get("jobName") == job_name