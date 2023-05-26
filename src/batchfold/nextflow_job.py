# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from attrs import define, field
from batchfold.batchfold_job import BatchFoldJob
from datetime import datetime
import logging

@define
class NextFlowJob(BatchFoldJob):
    """Define NextFlow Job"""

    assets_s3_uri: str = ""
    nextflow_script: str = ""
    params: dict = {}
    job_name: str = field(default="NextFlowJob" + datetime.now().strftime("%Y%m%d%s"))
    job_definition_name: str = field(default="NextflowJobDefinition")
    cpu: int = field(default=8)
    memory: int = field(default=15)
    gpu: int = field(default=0)

    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""
        command_list = [f"-i {self.assets_s3_uri}/:/home/"]
        command_list.extend(["nextflow", "run", self.nextflow_script])
        command_list.extend([f"--{key}=\'{value}\'" for key, value in self.params.items()])

        logging.info(f"Command is \n{command_list}")
        self.define_container_overrides(command_list, self.cpu, self.memory, self.gpu)
        return None