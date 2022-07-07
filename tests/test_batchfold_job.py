import pytest
import boto3
from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.batchfold_job import BatchFoldJob

@pytest.fixture()
def batch_environment():
    stack = BatchFoldEnvironment(boto_session = boto3.Session())
    assert "BatchEnvironment" in stack.stack_name
    return(stack)

def test_job_init():
    job = BatchFoldJob(
        job_name="MyNewJob",
        job_definition="MyJobDef",
        cpu = 6,
        memory = 9,
        gpu = 3,
        command = ["echo hello"],
        depends_on="MyOldJob"
    )

    assert job.job_name == "MyNewJob"
    assert job.job_definition == "MyJobDef"
    assert job.depends_on == "MyOldJob"
    assert job.container_overrides["command"] == ["echo hello"]
    assert len(job.container_overrides["resourceRequirements"]) == 3
    assert job.container_overrides["resourceRequirements"][0]["type"] == "VCPU"
    assert job.container_overrides["resourceRequirements"][0]["value"] == "6"
    assert job.container_overrides["resourceRequirements"][1]["type"] == "MEMORY"
    assert job.container_overrides["resourceRequirements"][1]["value"] == "9000"
    assert job.container_overrides["resourceRequirements"][2]["type"] == "GPU"
    assert job.container_overrides["resourceRequirements"][2]["value"] == "3"


def test_submit_batchfold_jobs(batch_environment):
    job_name = "TestJobNoOverrides"
    job_queue = batch_environment.queues["CPUOnDemandJobQueue"]
    job_definition = batch_environment.job_definitions["CPUFoldingJobDefinition"]

    new_job = BatchFoldJob(job_name=job_name, job_definition=job_definition)
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