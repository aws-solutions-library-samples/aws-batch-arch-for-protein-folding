from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.openfold_job import OpenFoldJob
import pytest
import boto3

@pytest.fixture()
def batch_environment():
    stack = BatchFoldEnvironment(boto_session = boto3.Session())
    assert "BatchEnvironment" in stack.stack_name
    return(stack)

def test_openfold_job_init():
    of_job = OpenFoldJob(
        job_definition = "Test",
        target_id = "T1084",
        fasta_s3_uri = "s3://aws-af-testing/T1084/fasta/T1084.fasta",
        output_s3_uri = "s3://aws-af-testing/T1084/outputs/",
        openfold_checkpoint_path = "/database/openfold_params/finetuning_ptm_1.pt",
    )

    assert of_job.target_id == "T1084"
    assert of_job.fasta_s3_uri == "s3://aws-af-testing/T1084/fasta/T1084.fasta"
    assert of_job.output_s3_uri == "s3://aws-af-testing/T1084/outputs/"
    assert of_job.template_mmcif_dir == "/database/pdb_mmcif/mmcif_files" 

def test_openfold_job_submission(batch_environment):

    job_name = "OpenFoldJob"
    job_queue = batch_environment.queues["G4dnJobQueue"]
    job_definition = batch_environment.job_definitions["OpenFoldJobDefinition"]

    of_job = OpenFoldJob(
        job_name = "OpenFoldJob",
        job_definition = job_definition,
        target_id = "T1084",
        fasta_s3_uri = "s3://aws-af-testing/T1084/fasta/T1084.fasta",
        output_s3_uri = "s3://aws-af-testing/T1084/outputs/",
        use_precomputed_alignments = "s3://aws-af-testing/T1084/msas/",
        config_preset = "finetuning_ptm",
        openfold_checkpoint_path = "/database/openfold_params/finetuning_ptm_1.pt",
        save_outputs = True        
    )

    response = job_queue.submit_job(of_job)
    assert job_name == response["jobName"] 
    
    job_description = of_job.describe_job()        
    assert job_name == job_description[0].get("jobName", [])

    pending_jobs = job_queue.list_jobs(
        valid_statuses = ["PENDING"]
    )
    assert bool(pending_jobs) is False

    job_list = job_queue.list_jobs()
    assert len(job_list) > 0

    job_info = [job for job in job_list if job.get("jobName", []) == job_name]
    assert job_info[0].get("jobDefinition") == job_definition