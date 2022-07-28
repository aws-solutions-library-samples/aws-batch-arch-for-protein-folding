from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.mmseqs2_job import MMseqs2Job
import pytest
import boto3
import time
from datetime import datetime

@pytest.fixture()
def batch_environment():
    stack = BatchFoldEnvironment(boto_session = boto3.Session())
    assert "BatchEnvironment" in stack.stack_name
    return(stack)

def test_mmseqs2_job_init():
    new_job = MMseqs2Job(
        target_id = "T1084",
        fasta_s3_uri = "s3://aws-af-testing/T1084/fasta/T1084.fasta",
        output_s3_uri = "s3://aws-af-testing/T1084/outputs/",
        mmseqs_database_path = "TESTA",
        pdb70_database_path = "TESTB"
    )

    assert new_job.job_definition_name == "MSAJobDefinition"
    assert new_job.target_id == "T1084"
    assert new_job.fasta_s3_uri == "s3://aws-af-testing/T1084/fasta/T1084.fasta"
    assert new_job.output_s3_uri == "s3://aws-af-testing/T1084/outputs/"
    assert new_job.mmseqs_database_path == "TESTA" 
    assert new_job.pdb70_database_path == "TESTB" 

def test_mmseqs2_job_submission(batch_environment):
    job_name = "MMseqs2Job" + datetime.now().strftime("%Y%m%d%s")
    job_queue_name = "GravitonSpotJobQueue"

    new_job = MMseqs2Job(
        job_name = job_name,
        target_id = "T1084",
        fasta_s3_uri = "s3://aws-af-testing/T1084/fasta/T1084.fasta",
        output_s3_uri = "s3://aws-af-testing/T1084/outputs/",
        cpu = 64,
        memory = 500
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