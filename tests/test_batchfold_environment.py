import pytest
import boto3
from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.batchfold_job import BatchFoldJob

@pytest.fixture()
def batch_environment():
    stack = BatchFoldEnvironment(boto_session = boto3.Session())
    assert "BatchEnvironment" in stack.stack_name
    return(stack)

def test_get_stack_outputs(batch_environment):
    outputs = batch_environment.get_stack_outputs()
    output_keys = list(outputs)
    output_keys.sort()
    assert output_keys == [
        "AlphaFold2JobDefinition",
        "CPUFoldingJobDefinition",
        "CPUOnDemandJobQueue",
        "CPUSpotJobQueue",
        "DownloadJobDefinition",
        "GPUJobQueue",
        "LaunchTemplate",
        "OpenFoldJobDefinition",
    ]

    assert "JobDefinition" in outputs["AlphaFold2JobDefinition"]
    assert "JobQueue" in outputs["CPUOnDemandJobQueue"]

def test_get_job_queue_names(batch_environment):
    outputs = batch_environment.get_stack_outputs(filter="JobQueue")
    output_keys = list(outputs)
    output_keys.sort()
    assert output_keys == [
        "CPUOnDemandJobQueue",
        "CPUSpotJobQueue",
        "GPUJobQueue",
    ]

def test_get_job_definition_names(batch_environment):
    outputs = batch_environment.get_stack_outputs(filter="JobDefinition")
    output_keys = list(outputs)
    output_keys.sort()
    job_defs = list(batch_environment.job_definitions)
    job_defs.sort()
    assert output_keys == job_defs == [
        "AlphaFold2JobDefinition",
        "CPUFoldingJobDefinition",
        "DownloadJobDefinition",
        "OpenFoldJobDefinition",
    ]

def test_get_job_queue_objects(batch_environment):
    assert len(batch_environment.queues) == 3
    assert batch_environment.queues["CPUOnDemandJobQueue"].type == "CPUOnDemandJobQueue"
    assert batch_environment.queues["CPUOnDemandJobQueue"].jobs == []


def test_submit_batchfold_jobs(batch_environment):
    job_name = "TestJobNoOverrides"
    job_queue = batch_environment.queues["CPUOnDemandJobQueue"]
    job_definition = batch_environment.job_definitions["CPUFoldingJobDefinition"]

    new_job = BatchFoldJob(job_name, job_definition)
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