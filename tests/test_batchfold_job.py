import pytest
import boto3
from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.batchfold_job import BatchFoldJob


@pytest.fixture()
def batch_environment():
    stack = BatchFoldEnvironment(boto_session=boto3.Session())
    assert "BatchEnvironment" in stack.stack_name
    return stack


def test_job_init():
    job = BatchFoldJob(
        job_name="MyNewJob",
        job_definition_name="MyJobDef",
        cpu=6,
        memory=9,
        gpu=3,
        command=["echo hello"],
    )

    assert job.job_name == "MyNewJob"
    assert job.job_definition_name == "MyJobDef"
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
    job_queue_name = "CPUOnDemandJobQueue"

    new_job = BatchFoldJob(job_name=job_name)
    batch_environment.submit_job(job=new_job, job_queue_name=job_queue_name)
    assert job_name == batch_environment.last_submission.job_name

    job_description = new_job.describe_job()
    assert job_name == job_description[0].get("jobName", [])

    pending_jobs = batch_environment.list_jobs(valid_statuses=["PENDING"])
    assert bool(pending_jobs[job_queue_name]) is False

    job_dict = batch_environment.list_jobs()
    job_list = job_dict[job_queue_name]
    assert len(job_list) > 0

    job_info = [job for job in job_list if job.get("jobName", []) == job_name]
    assert (
        job_info[0].get("jobDefinition")
        == batch_environment.job_definitions["CPUFoldingJobDefinition"]
    )
    assert job_info[0].get("jobName") == job_name


def test_submit_chained_jobs(batch_environment):

    job_queue_name = "CPUOnDemandJobQueue"

    job_1 = BatchFoldJob(job_name="ChainedJob1")
    job_2 = BatchFoldJob(job_name="ChainedJob2")
    batch_environment \
        .submit_job(job=job_1, job_queue_name=job_queue_name) \
        .submit_job(job=job_2, job_queue_name=job_queue_name)

    assert job_1.describe_job()[0].get("jobID", [])  != ""
    assert job_2.describe_job()[0].get("jobID", [])  != ""

    job_3 = BatchFoldJob(job_name="ChainedJob3")
    job_4 = BatchFoldJob(job_name="ChainedJob4")
    batch_environment \
        .submit_job(job=job_3, job_queue_name=job_queue_name) \
        .submit_job(job=job_4, job_queue_name=job_queue_name, dependent=True)

    assert job_3.describe_job()[0].get("jobID", [])  != ""
    assert job_4.describe_job()[0].get("jobID", [])  != ""
