from attrs import define, field
from typing import List, Dict
import boto3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

@define
class BatchFoldJob:
    """Base class to submit protein analysis job to AWS"""

    boto_session: boto3.session.Session = field(
        default=boto3.DEFAULT_SESSION or boto3.Session()
    )
    id: str = field(init=False)
    container_overrides: Dict = field(init=False)  # Depends on cpu
    def define_container_overrides(self, command, cpu=4, memory=16, gpu=0) -> None:
        """Convert the job parameters into a container_overrides object for AWS Batch"""
        self.container_overrides = {
            "command": command,
            "resourceRequirements": [
                {"value": str(cpu), "type": "VCPU"},
                {"value": str(memory * 1000), "type": "MEMORY"},
            ],
        }
        if gpu > 0:
            self.container_overrides["resourceRequirements"].append(
                {"value": str(gpu), "type": "GPU"}
            )
        return None

    def describe_job(self) -> List:
        """Describe AWS Batch jobs"""

        batch = self.boto_session.client("batch")
        return batch.describe_jobs(jobs=[self.id]).get("jobs", [])
