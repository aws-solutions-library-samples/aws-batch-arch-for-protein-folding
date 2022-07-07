import boto3
from attrs import define, field
from typing import List, Dict
from batchfold.batchfold_job_queue import JobQueue

@define
class BatchFoldEnvironment:
    """Manage protein analysis environments on AWS Batch.

    This class provides methods for submitting jobs to AWS Batch.

    Attributes:
        boto_session (boto3.session.Session): The underlying Boto3 session which AWS service
            calls are delegated to (default: None). If not provided, one is created with
            default AWS configuration chain.
        stack_name (str): Name of the CloudFormation stack to use. If not provided, the
            most recent stack will be used.
        queues (dict): Dictionary of JobQueue objects used to submit and track analysis jobs.
        job_definitions (list): List of valid job definitions for the environment.

    """
    
    boto_session: boto3.session.Session = boto3.DEFAULT_SESSION or boto3.Session()
    stack_name: str = field(kw_only=True)
    queues: dict = field(kw_only=True)
    job_definitions: dict = field(kw_only=True)

    @stack_name.default
    def _get_latest_stack(self) -> str:
        """Get the latest batchfold cloudformation stack."""

        cfn = self.boto_session.client("cloudformation")
        batchfold_stacks = []
        stack_status_filter=["CREATE_COMPLETE", "ROLLBACK_COMPLETE", "UPDATE_COMPLETE"]
        for stack in cfn.list_stacks(StackStatusFilter=stack_status_filter)[
            "StackSummaries"
        ]:
            if "batch-protein-folding-cfn-batch.yaml" in stack.get(
                "TemplateDescription", []
            ):
                batchfold_stacks.append(stack)
        return batchfold_stacks[0].get("StackName", [])

    @queues.default
    def _get_queue_dict(self) -> Dict:
        """Create new JobQueue instances and add them to the BatchFold environment instance."""
        
        queues = self.get_stack_outputs(filter="JobQueue")
        queue_dict = {}
        for queue_type, queue_name in queues.items():
            queue = JobQueue(queue_type, queue_name)
            queue_dict[queue_type] = queue
        return queue_dict

    @job_definitions.default
    def _get_job_definitions_names(self) -> Dict:
        """Get the valid job definition names and add them to the BatchFold Environment instance."""
        
        job_definitions = self.get_stack_outputs(filter="JobDefinition")
        return job_definitions

    def get_stack_outputs(self, filter: str = "") -> Dict:
        """Get a dict of the cloudformation stack outputs, optionally with key filtered by a string"""

        cfn = self.boto_session.client("cloudformation")
        output_list = cfn.describe_stacks(StackName=self.stack_name).get("Stacks")[0].get("Outputs")
        output_dict = {output["OutputKey"]:output["OutputValue"] for output in output_list if filter in output["OutputKey"]}
        return output_dict
    
