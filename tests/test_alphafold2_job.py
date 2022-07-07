from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.alphafold2_job import AlphaFold2Job
import pytest
import boto3
from datetime import datetime

@pytest.fixture()
def batch_environment():
    stack = BatchFoldEnvironment(boto_session = boto3.Session())
    assert "BatchEnvironment" in stack.stack_name
    return(stack)

def test_alphafold_2_job_init():
    new_job = AlphaFold2Job(
        job_definition = "Test",
        target_id = "T1084",
        fasta_s3_uri = "s3://aws-af-testing/T1084/fasta/T1084.fasta",
        output_s3_uri = "s3://aws-af-testing/T1084/outputs/",
        max_template_date = "2022-05-31",
    )

    assert new_job.target_id == "T1084"
    assert new_job.fasta_s3_uri == "s3://aws-af-testing/T1084/fasta/T1084.fasta"
    assert new_job.output_s3_uri == "s3://aws-af-testing/T1084/outputs/"
    assert new_job.max_template_date == "2022-05-31" 

def test_alphafold_2_job_submission(batch_environment):

    job_name = "AlphaFoldJ2ob" + datetime.now().strftime("%Y%m%d%s")
    job_queue = batch_environment.queues["G4dnJobQueue"]
    job_definition = batch_environment.job_definitions["AlphaFold2JobDefinition"]

    new_job = AlphaFold2Job(
        job_name = job_name,
        job_definition = job_definition,
        target_id = "T1084",
        fasta_s3_uri = "s3://aws-af-testing/T1084/fasta/T1084.fasta",
        msa_s3_uri = "s3://aws-af-testing/T1084/msas/",
        output_s3_uri = "s3://aws-af-testing/T1084/outputs/",
        use_precomputed_msas = True,
        model_preset = "monomer"
    )

    response = job_queue.submit_job(new_job)
    assert job_name == response["jobName"] 
    
    job_description = new_job.describe_job()        
    assert job_name == job_description[0].get("jobName", [])

    pending_jobs = job_queue.list_jobs(
        valid_statuses = ["PENDING"]
    )
    assert bool(pending_jobs) is False

    job_list = job_queue.list_jobs()
    assert len(job_list) > 0

    job_info = [job for job in job_list if job.get("jobName", []) == job_name]
    assert job_info[0].get("jobDefinition") == job_definition