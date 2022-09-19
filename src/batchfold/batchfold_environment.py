from __future__ import annotations

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

import boto3
from attrs import define, field
from typing import List, Dict
from batchfold.batchfold_job_queue import JobQueue, JobSubmission
from batchfold.batchfold_job import BatchFoldJob


@define
class BatchFoldEnvironment:
    """Manage protein analysis environments on AWS Batch.

    This class provides methods for submitting jobs to AWS Batch.

    Attributes:
        boto_session (boto3.session.Session): The underlying Boto3 session which AWS service
            calls are delegated to (default: None). If not provided, one is created with
            default AWS configuration chain.
        stack_id (str): ID of the CloudFormation stack to use. If not provided, the
            most recent stack will be used.
        queues (dict): Dictionary of JobQueue objects used to submit and track analysis jobs.
        job_definitions (list): List of valid job definitions for the environment.

    """

    boto_session: boto3.session.Session = boto3.DEFAULT_SESSION or boto3.Session()
    stack_id: str = field(kw_only=True)
    nested_stacks: List = field(kw_only=True)
    stack_outputs: List = field(kw_only=True)
    queues: dict = field(kw_only=True)
    job_definitions: dict = field(kw_only=True)
    default_bucket: str = field(kw_only=True)

    @stack_id.default
    def load_latest_root_stack(self) -> str:
        """Get the latest batchfold cloudformation stack."""

        cfn = self.boto_session.client("cloudformation")
        batchfold_stacks = []
        stack_status_filter = [
            "CREATE_COMPLETE",
            "ROLLBACK_COMPLETE",
            "UPDATE_COMPLETE",
        ]
        for stack in cfn.list_stacks(StackStatusFilter=stack_status_filter)["StackSummaries"]:
            if "batch-protein-folding-cfn-root.yaml" in stack.get("TemplateDescription", []):
                batchfold_stacks.append((stack["CreationTime"], stack["StackId"]))
        return sorted(batchfold_stacks, key = lambda x: x[0], reverse=True)[0][1]

    @nested_stacks.default
    def load_nested_stacks(self) -> List:
        cfn = self.boto_session.client("cloudformation")
        resources = cfn.list_stack_resources(StackName=self.stack_id).get("StackResourceSummaries", [])
        return [ x.get("PhysicalResourceId", []) for x in resources if x.get("ResourceType", []) ==  "AWS::CloudFormation::Stack" ]

    @stack_outputs.default
    def load_stack_outputs(self) -> List:
        cfn = self.boto_session.client("cloudformation")
        stack_outputs = []
        for nested_stack in self.nested_stacks:
            stack_outputs.extend(cfn.describe_stacks(StackName=nested_stack).get("Stacks")[0].get("Outputs", []))
        return(stack_outputs)

    def get_stack_outputs(self, filter: str = "") -> Dict:
        """Get a dict of the cloudformation stack outputs, optionally with key filtered by a string"""
        output_dict = {
            output["OutputKey"]: output["OutputValue"]
            for output in self.stack_outputs
            if filter in output["OutputKey"]
        }
        return output_dict

    @queues.default
    def load_queues(self) -> Dict:
        """Create new JobQueue instances and add them to the BatchFold environment instance."""

        queues = self.get_stack_outputs(filter="JobQueue")
        queue_dict = {}
        for queue_name, queue_id in queues.items():
            queue = JobQueue(queue_name, queue_id, boto_session=self.boto_session)
            queue_dict[queue_name] = queue
        return queue_dict

    def list_job_queue_names(self) -> List:
        """List the names of all available job queues."""
        queue_names = list(self.queues.keys())
        queue_names.sort()
        return queue_names

    @job_definitions.default
    def load_job_definitions(self) -> Dict:
        """Get the valid job definition names and add them to the BatchFold Environment instance."""

        job_definitions = self.get_stack_outputs(filter="JobDefinition")
        return job_definitions

    def list_job_definition_names(self) -> List:
        """List the names of all available job definitions."""
        names = list(self.job_definitions.keys())
        names.sort()
        return names

    @default_bucket.default
    def load_default_bucket(self):
        return self.get_stack_outputs().get("S3BucketName", [])

    def submit_job(
        self,
        job: BatchFoldJob,
        job_queue_name: str,
        depends_on: List[JobSubmission] = None,
    ) -> BatchFoldEnvironment:
        """Submit job to specified job queue."""

        job_queue = self.queues[job_queue_name]
        job_definition = self.job_definitions[job.job_definition_name]
        return job_queue.submit_job(job, job_definition, depends_on)

    def list_jobs(
        self,
        statuses: List = [
            "SUBMITTED",
            "PENDING",
            "RUNNABLE",
            "STARTING",
            "RUNNING",
            "SUCCEEDED",
            "FAILED",
        ],
        queues: List = None,
    ) -> Dict:
        """List jobs on all job queues."""
        queues = queues or self.list_job_queue_names()
        all_jobs = {}
        for queue in [queue for (name, queue) in self.queues.items() if name in queues]:
            jobs = queue.list_jobs(statuses=statuses)
            all_jobs[queue.name] = jobs

        return all_jobs

