import boto3
from attrs import define
from typing import List, Dict
import time
from batchfold.batchfold_job import BatchFoldJob

@define
class JobQueue:

    type: str
    name: str
    boto_session: boto3.session.Session = boto3.DEFAULT_SESSION or boto3.Session()
    jobs: list = []

    def list_jobs(
        self, 
        valid_statuses: List = ['SUBMITTED','PENDING', 'RUNNABLE', 'STARTING', 'RUNNING', 'SUCCEEDED', 'FAILED']) -> List:
        """List jobs on the Job Queue instance."""
        
        batch = self.boto_session.client("batch")
        jobs = batch.list_jobs(jobQueue = self.name, filters = [{"name":"BEFORE_CREATED_AT", "values":[str(int(time.time()))]}])
        filtered_jobs = [job for job in jobs.get("jobSummaryList", []) if job.get("status") in valid_statuses]
        return filtered_jobs

    def submit_job(self, 
        job: BatchFoldJob
        ) -> Dict:
        """Submit Batch job"""
        batch = self.boto_session.client("batch")

        if job.container_overrides and job.depends_on:
            response = batch.submit_job(
                jobName=job.name,
                jobQueue=self.name,
                jobDefinition=job.definition,
                containerOverrides=job.container_overrides,
                dependsOn=job.depends_on
            )
        elif job.container_overrides and not job.depends_on:
            response = batch.submit_job(
                jobName=job.name,
                jobQueue=self.name,
                jobDefinition=job.definition,
                containerOverrides=job.container_overrides
            )
        elif job.depends_on and not job.container_overrides:
            response = batch.submit_job(
                jobName=job.name,
                jobQueue=self.name,
                jobDefinition=job.definition,
                dependsOn=job.depends_on
            )     
        else:
            response = batch.submit_job(
                jobName=job.name,
                jobQueue=self.name,
                jobDefinition=job.definition
            )
        job.id = response.get("jobId", [])
        return response