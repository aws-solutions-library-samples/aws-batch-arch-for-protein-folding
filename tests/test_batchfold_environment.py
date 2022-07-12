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
        "G4dnJobQueue",
        "GravitonJobQueueSpot",
        "LaunchTemplate",
        "OpenFoldJobDefinition",
    ]

    assert "JobDefinition" in outputs["AlphaFold2JobDefinition"]
    assert "JobQueue" in outputs["CPUOnDemandJobQueue"]
    assert "JobQueue" in outputs["G4dnJobQueue"]

def test_get_job_queue_names(batch_environment):
    outputs = batch_environment.get_stack_outputs(filter="JobQueue")
    output_keys = list(outputs)
    output_keys.sort()
    assert output_keys == [
        "CPUOnDemandJobQueue",
        "CPUSpotJobQueue",
        "G4dnJobQueue",
        "GravitonJobQueueSpot"
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
    assert len(batch_environment.queues) == 4
    assert batch_environment.queues["CPUOnDemandJobQueue"].type == "CPUOnDemandJobQueue"
    assert batch_environment.queues["CPUOnDemandJobQueue"].jobs == []

