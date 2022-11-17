# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.rfdesign_job import RFDesignHallucinateJob
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

def test_rfdesign_job_init():
    bucket = os.getenv("TEST_BUCKET")
    job_name = "RFDesignJob" + datetime.now().strftime("%Y%m%d%s")
    params = {
        "contigs": "A163-181",
        "len": 60,
        "save_pdb": True
    }
    new_job = RFDesignHallucinateJob(
        boto_session=boto_session,
        job_name = job_name,
        memory=61,
        target_id = "5TPN",
        input_s3_uri = f"s3://{bucket}/rfdesign/",
        output_s3_uri = f"s3://{bucket}/rfdesign/5tpn/outputs/",
        pdb = "rsvf-v_5tpn.pdb",
        params = params
    )
    assert new_job.job_definition_name == "RFDesignJobDefinition"
    assert new_job.target_id == "5TPN"
    assert new_job.input_s3_uri == f"s3://{bucket}/rfdesign/"
    assert new_job.output_s3_uri == f"s3://{bucket}/rfdesign/5tpn/outputs/"
    assert new_job.pdb == "rsvf-v_5tpn.pdb"
    assert new_job.params['contigs'] == "A163-181"
    assert new_job.params['len'] == 60
    assert new_job.params['save_pdb'] == True

def test_rfdesign_job_submission(batch_environment):

    job_name = "RFDesignJob" + datetime.now().strftime("%Y%m%d%s")
    job_queue_name = "G4dnJobQueue"
    bucket = os.getenv("TEST_BUCKET")
    params = {
        "contigs": "A163-181",
        "len": 60,
        "steps": "g3",
        "w_rep": 1,
        "rep_sigma": 4,
        "rep_pdb": "rsvf-v_5tpn_receptor_frag.pdb",
        "w_atr": 10,
        "atr_sigma": 6,
        "atr_pdb": "rsvf-v_5tpn_receptor_frag.pdb",
        "save_pdb": True,
        "n_bkg": 5
    }
    new_job = RFDesignHallucinateJob(
        boto_session=boto_session,
        job_name = job_name,
        memory=61,
        target_id = "5TPN",
        input_s3_uri = f"s3://{bucket}/rfdesign/",
        output_s3_uri = f"s3://{bucket}/rfdesign/5tpn/outputs/",
        pdb = "rsvf-v_5tpn.pdb",
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