from attrs import define, field
from typing import List, Dict
import boto3
from datetime import datetime

@define
class BatchFoldJob:
    """ Submit protein analysis job to AWS """

    job_name: str = datetime.now().strftime("%Y%m%dT%H%M%S")
    job_definition: str = field(kw_only=True)
    cpu: int = 4
    memory: int = 16
    gpu: int = 0
    command: List = ["echo hello"]
    depends_on: str = ""    
    container_overrides: Dict = field(kw_only=True) 
    boto_session: boto3.session.Session = boto3.DEFAULT_SESSION or boto3.Session()
    id: str = field(init=False)

    @container_overrides.default
    def _define_container_overrides(self):
        """Convert the job parameters into a container_overrides object for AWS Batch"""
        container_overrides = {
            "command": self.command,
            "resourceRequirements": [
                {"value": str(self.cpu), "type": "VCPU"},
                {"value": str(self.memory * 1000), "type": "MEMORY"},
            ]
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
