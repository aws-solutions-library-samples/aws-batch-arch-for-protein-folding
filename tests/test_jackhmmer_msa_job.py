from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.jackhmmer_job import JackhmmerJob
import pytest
import boto3
from datetime import datetime


@pytest.fixture()
def batch_environment():
    stack = BatchFoldEnvironment(boto_session = boto3.Session())
    assert "BatchEnvironment" in stack.stack_name
    return(stack)

def test_jackhmmer_job_init():
    new_job = JackhmmerJob(
        job_definition="Test",
        target_id = "T1084",
        fasta_s3_uri = "s3://aws-af-testing/T1084/fasta/T1084.fasta",
        output_s3_uri = "s3://aws-af-testing/T1084/outputs/",
        bfd_database_path = "TESTA",
        uniclust30_database_path = "TESTB"
    )

    assert new_job.target_id == "T1084"
    assert new_job.fasta_s3_uri == "s3://aws-af-testing/T1084/fasta/T1084.fasta"
    assert new_job.output_s3_uri == "s3://aws-af-testing/T1084/outputs/"
    assert new_job.bfd_database_path == "TESTA" 
    assert new_job.uniclust30_database_path == "TESTB" 

def test_jackhmmer_job_submission(batch_environment):

    job_name = "JackhmmerJob" + datetime.now().strftime("%Y%m%d%s")
    job_queue = batch_environment.queues["GravitonJobQueueSpot"]
    job_definition = batch_environment.job_definitions["MSAJobDefinition"]

    new_job = JackhmmerJob(
        job_name = job_name,
        job_definition = job_definition,
        target_id = "T1084",
        fasta_s3_uri = "s3://aws-af-testing/T1084/fasta/T1084.fasta",
        output_s3_uri = "s3://aws-af-testing/T1084/outputs/",
        bfd_database_path = "small_bfd/tiny.fasta",
        mgnify_database_path = None,
        pdb70_database_path = None,
        uniprot_database_path = None,
        uniref90_database_path = None,
        use_small_bfd = True,
        cpus = 4
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