# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

import pytest
import boto3
from batchfold.batchfold_environment import BatchFoldEnvironment

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
        "G4dnJobQueue",
        "GravitonOnDemandJobQueue",
        "GravitonSpotJobQueue",
        "LaunchTemplate",
        "MSAJobDefinition",
        "OpenFoldJobDefinition",
        "S3BucketName"
    ]

    assert "JobDefinition" in outputs["AlphaFold2JobDefinition"]
    assert "JobQueue" in outputs["CPUOnDemandJobQueue"]
    assert "JobQueue" in outputs["G4dnJobQueue"]
    assert "batchfolds3bucket" in outputs["S3BucketName"]

def test_get_job_queue_names(batch_environment):
    job_queue_names = batch_environment.list_job_queue_names()
    assert job_queue_names == [
        "CPUOnDemandJobQueue",
        "CPUSpotJobQueue",
        "G4dnJobQueue",
        "GravitonOnDemandJobQueue",
        "GravitonSpotJobQueue"
    ]

def test_get_job_definition_names(batch_environment):
    job_def_names = batch_environment.list_job_definition_names()
    assert job_def_names == [
        "AlphaFold2JobDefinition",
        "CPUFoldingJobDefinition",
        "MSAJobDefinition",
        "OpenFoldJobDefinition",
    ]

def test_get_job_queue_objects(batch_environment):
    assert len(batch_environment.queues) == 5
    assert batch_environment.queues["CPUOnDemandJobQueue"].name == "CPUOnDemandJobQueue"

def test_get_default_bucket(batch_environment):
    assert "batchfolds3bucket" in batch_environment.default_bucket