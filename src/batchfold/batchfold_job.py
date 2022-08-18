from attrs import define, field
from typing import List, Dict
import boto3
from datetime import datetime

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0


@define
class BatchFoldJob:
    """Base class to submit protein analysis job to AWS"""

    boto_session: boto3.session.Session = field(
        kw_only=True, default=boto3.DEFAULT_SESSION or boto3.Session()
    )
    command: List = field(kw_only=True, default=["echo hello"])
    cpu: int = field(kw_only=True, default=4)
    job_name: str = field(
        kw_only=True, default=datetime.now().strftime("%Y%m%dT%H%M%S")
    )
    job_definition_name: str = field(kw_only=True, default="CPUFoldingJobDefinition")
    gpu: int = field(kw_only=True, default=0)
    id: str = field(init=False)
    memory: int = field(kw_only=True, default=16)
    container_overrides: Dict = field(kw_only=True)  # Depends on cpu

    @container_overrides.default
    def _define_container_overrides(self):
        """Convert the job parameters into a container_overrides object for AWS Batch"""
        container_overrides = {
            "command": self.command,
            "resourceRequirements": [
                {"value": str(self.cpu), "type": "VCPU"},
                {"value": str(self.memory * 1000), "type": "MEMORY"},
            ],
        }

        if self.gpu > 0:
            container_overrides["resourceRequirements"].append(
                {"value": str(self.gpu), "type": "GPU"}
            )
        return container_overrides

    def describe_job(self) -> List:
        """Describe AWS Batch jobs"""

        batch = self.boto_session.client("batch")
        return batch.describe_jobs(jobs=[self.id]).get("jobs", [])
