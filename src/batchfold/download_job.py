# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from attrs import define
from batchfold.batchfold_job import BatchFoldJob
import logging


@define
class DownloadJob(BatchFoldJob):
    """Define Download MSA Creation Job"""

    job_definition_name: str = "DownloadJobDefinition"
    script: str = "./scripts/download_test.sh"
    data_dir: str = "/database"

    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""

        command_string = f"{self.script} {self.data_dir}"

        logging.info(f"Command is \n{command_string}")
        self.container_overrides["command"] = [command_string]

        return None
