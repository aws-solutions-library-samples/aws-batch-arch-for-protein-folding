# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from attrs import define, field
from batchfold.batchfold_job import BatchFoldJob
from datetime import datetime
import logging

@define
class DownloadJob(BatchFoldJob):
    """Define Download MSA Creation Job"""

    script: str = "./scripts/download_test.sh"
    data_dir: str = "/database"
    job_name: str = field(default="DownloadJob" + datetime.now().strftime("%Y%m%d%s"))
    job_definition_name: str = field(default="DownloadJobDefinition")
    cpu: int = field(default=4)
    memory: int = field(default=16)
    gpu: int = field(default=0)

    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""

        command_string = f"{self.script} {self.data_dir}"
        logging.info(f"Command is \n{command_string}")
        self.define_container_overrides([command_string], self.cpu, self.memory, self.gpu)
        return None