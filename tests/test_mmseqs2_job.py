# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.mmseqs2_job import MMSeqs2Job
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

def test_mmseqs2_job_init():
    bucket = os.getenv("TEST_BUCKET")
    new_job = MMSeqs2Job(
        boto_session = boto_session,
        target_id = "T1084",
        fasta_s3_uri = f"s3://{bucket}/T1084/fasta/T1084.fasta",
        output_s3_uri = f"s3://{bucket}/T1084/outputs/",
        databases = ["uniref90", "swissprot"],
        databases_paths_overwrite = {
            "uniref90": "mmseqs2_dbs/uniref90/uniref90",
            "swissprot": "mmseqs2_dbs/swissprot/swissprot",
        },
    )
    assert new_job.job_definition_name == "MMSeqs2JobDefinition"
    assert new_job.target_id == "T1084"
    assert new_job.fasta_s3_uri == f"s3://{bucket}/T1084/fasta/T1084.fasta"
    assert new_job.output_s3_uri == f"s3://{bucket}/T1084/outputs/"
    assert new_job.databases == ["uniref90", "swissprot"]
    assert new_job.databases_paths_overwrite["uniref90"] == "mmseqs2_dbs/uniref90/uniref90"

def test_mmseqs2_gpu_job_init():
    bucket = os.getenv("TEST_BUCKET")
    new_job = MMSeqs2Job(
        boto_session = boto_session,
        target_id = "T1084",
        fasta_s3_uri = f"s3://{bucket}/T1084/fasta/T1084.fasta",
        output_s3_uri = f"s3://{bucket}/T1084/outputs/",
        gpu = 1,
    )
    assert new_job.job_definition_name == "MMSeqs2GPUJobDefinition"
    assert new_job.gpu == 1

def test_mmseqs2_command_construction():
    bucket = os.getenv("TEST_BUCKET")
    new_job = MMSeqs2Job(
        boto_session = boto_session,
        target_id = "T1084",
        fasta_s3_uri = f"s3://{bucket}/T1084/fasta/T1084.fasta",
        output_s3_uri = f"s3://{bucket}/T1084/outputs/",
        databases = ["uniref90"],
        evalue = 0.001,
        cov = 0.8,
        cov_mode = 1,
        max_seqs = 5000,
        sensitivity = 9.0,
        num_iterations = 5,
        output_format = "sto",
        write_hits_table = False,
        cpu = 8,
    )
    command_str = " ".join(new_job.container_overrides["command"])
    assert f"-i s3://{bucket}/T1084/fasta/T1084.fasta:" in command_str
    assert "run_mmseqs_msa.sh" in command_str
    assert "--evalue 0.001" in command_str
    assert "--cov 0.8" in command_str
    assert "--cov_mode 1" in command_str
    assert "--max_seqs 5000" in command_str
    assert "--s 9.0" in command_str
    assert "--num_iterations 5" in command_str
    assert "--output_format sto" in command_str
    assert "--write_hits_table 0" in command_str
    assert "--threads 8" in command_str
    assert "uniref90:uniref90/uniref90" in command_str

def test_mmseqs2_database_override():
    bucket = os.getenv("TEST_BUCKET")
    new_job = MMSeqs2Job(
        boto_session = boto_session,
        target_id = "T1084",
        fasta_s3_uri = f"s3://{bucket}/T1084/fasta/T1084.fasta",
        output_s3_uri = f"s3://{bucket}/T1084/outputs/",
        databases = ["uniref90", "mgnify"],
        databases_paths_overwrite = {
            "uniref90": "custom/uniref90/path",
            "mgnify": "custom/mgnify/path",
        },
    )
    command_str = " ".join(new_job.container_overrides["command"])
    assert "uniref90:custom/uniref90/path" in command_str
    assert "mgnify:custom/mgnify/path" in command_str

def test_mmseqs2_resource_requirements():
    bucket = os.getenv("TEST_BUCKET")
    new_job = MMSeqs2Job(
        boto_session = boto_session,
        target_id = "T1084",
        fasta_s3_uri = f"s3://{bucket}/T1084/fasta/T1084.fasta",
        output_s3_uri = f"s3://{bucket}/T1084/outputs/",
    )
    resources = new_job.container_overrides["resourceRequirements"]
    vcpu = next(r for r in resources if r["type"] == "VCPU")
    memory = next(r for r in resources if r["type"] == "MEMORY")
    assert vcpu["value"] == "16"
    assert memory["value"] == "31000"

def test_mmseqs2_gpu_resource_requirements():
    bucket = os.getenv("TEST_BUCKET")
    new_job = MMSeqs2Job(
        boto_session = boto_session,
        target_id = "T1084",
        fasta_s3_uri = f"s3://{bucket}/T1084/fasta/T1084.fasta",
        output_s3_uri = f"s3://{bucket}/T1084/outputs/",
        gpu = 1,
    )
    resources = new_job.container_overrides["resourceRequirements"]
    gpu = next((r for r in resources if r["type"] == "GPU"), None)
    assert gpu is not None
    assert gpu["value"] == "1"

def test_mmseqs2_job_submission(batch_environment):
    job_name = "MMSeqs2Job" + datetime.now().strftime("%Y%m%d%s")
    job_queue_name = "CPUOnDemandJobQueue"
    bucket = os.getenv("TEST_BUCKET")
    new_job = MMSeqs2Job(
        boto_session = boto_session,
        job_name = job_name,
        target_id = "T1084",
        fasta_s3_uri = f"s3://{bucket}/T1084/fasta/T1084.fasta",
        output_s3_uri = f"s3://{bucket}/T1084/outputs/",
        cpu = 4,
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
    assert job_info[0].get("jobDefinition") == batch_environment.job_definitions["MMSeqs2JobDefinition"]
    assert job_info[0].get("jobName") == job_name
