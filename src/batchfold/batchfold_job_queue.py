from __future__ import annotations

import boto3
from attrs import define
from typing import List, Dict
import time
from batchfold.batchfold_job import BatchFoldJob

@define
class JobQueue:

    name: str
    id: str
    boto_session: boto3.session.Session = boto3.DEFAULT_SESSION or boto3.Session()
    jobs: list = []

    def list_jobs(
        self, 
        statuses: List = ['SUBMITTED','PENDING', 'RUNNABLE', 'STARTING', 'RUNNING', 'SUCCEEDED', 'FAILED']) -> List:
        """List jobs on the Job Queue instance."""
        
        batch = self.boto_session.client("batch")
        jobs = batch.list_jobs(jobQueue = self.id, filters = [{"name":"BEFORE_CREATED_AT", "values":[str(int(time.time()))]}])
        filtered_jobs = [job for job in jobs.get("jobSummaryList", []) if job.get("status") in statuses]
        return filtered_jobs

    def submit_job(self, 
        job: BatchFoldJob,
        job_definition: str,
        depends_on: List[JobSubmission]
        ) -> Dict:
        """Submit Batch job"""
        batch = self.boto_session.client("batch")

        if job.container_overrides and depends_on:
            response = batch.submit_job(
                jobName=job.job_name,
                jobQueue=self.id,
                jobDefinition=job_definition,
                containerOverrides=job.container_overrides,
                dependsOn=[{"jobId":job.job_id, "type":"Sequential"} for job in depends_on]
            )
        elif job.container_overrides and not depends_on:
            response = batch.submit_job(
                jobName=job.job_name,
                jobQueue=self.id,
                jobDefinition=job_definition,
                containerOverrides=job.container_overrides
            )
        elif depends_on and not job.container_overrides:
            response = batch.submit_job(
                jobName=job.job_name,
                jobQueue=self.id,
                jobDefinition=job_definition,
                # dependsOn=[{"jobId": depends_on.job_id, "type": "Sequential"}]
                dependsOn=[{"jobId":job.job_id, "type":"Sequential"} for job in depends_on]                
            )     
        else:
            response = batch.submit_job(
                jobName=job.job_name,
                jobQueue=self.id,
                jobDefinition=job_definition
            )
        job.id = response.get("jobId", [])

        submission = JobSubmission(
            job_arn = response.get("jobArn", []),
            job_name = response.get("jobName", []),
            job_id = response.get("jobId", []),
            job_queue = self.id,
            job_definition = job_definition,
            job_type = type(job).__name__,
            depends_on = depends_on
        )

        return submission

@define
class JobSubmission:
    job_arn: str
    job_name: str
    job_id: str
    job_queue: str
    job_definition: str
    job_type: str
    depends_on: List