# import pytest
# import os
# import boto3
# from datetime import datetime
# from taskcat.testing import CFNTest
# from batchfold.batchfold_environment import BatchFoldEnvironment
# from batchfold.batchfold_job import BatchFoldJob

# @pytest.mark.skip(reason="Under development")
# def test_taskcat():
#     test = CFNTest.from_file()
#     with test as stacks:
#         # Calling 'with' or 'test.run()' will deploy the stacks.
#         for stack in stacks:

#             print(f"Testing {stack.name}")
#             assert "AlphaFold2JobDefinition" in stack.outputs["AlphaFold2JobDefinition"]
#             assert "CPUFoldingJobDefinition" in stack.outputs["CPUFoldingJobDefinition"]
#             assert "CPUOnDemandJobQueue" in stack.outputs["CPUOnDemandJobQueue"]
#             assert "CPUSpotJobQueue" in stack.outputs["CPUSpotJobQueue"]
#             assert "DownloadJobDefinition" in stack.outputs["DownloadJobDefinition"]
#             assert "G4dnJobQueue" in stack.outputs["G4dnJobQueue"]
#             assert "GravitonSpotJobQueue" in stack.outputs["GravitonSpotJobQueue"]
#             assert "LaunchTemplate" in stack.outputs["LaunchTemplate"]
#             assert "MSAJobDefinition" in stack.outputs["MSAJobDefinition"]
#             assert "OpenFoldJobDefinition" in stack.outputs["OpenFoldJobDefinition"]
#             assert "S3BucketName" in stack.outputs["S3BucketName"]

#             batch_environment = BatchFoldEnvironment(boto_session=boto3.Session())
#             job_name = "BatchFoldTestJobNoOverrides" + datetime.now().strftime("%Y%m%d%s")
#             job_queue_name = "CPUOnDemandJobQueue"
        
#             new_job = BatchFoldJob(job_name=job_name)
#             submission = batch_environment.submit_job(job=new_job, job_queue_name=job_queue_name)
#             assert job_name == submission.job_name

#             job_description = new_job.describe_job()
#             assert job_name == job_description[0].get("jobName", [])