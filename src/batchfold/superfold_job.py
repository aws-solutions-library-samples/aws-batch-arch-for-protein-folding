# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from attrs import define, field
from batchfold.batchfold_job import BatchFoldJob
from datetime import datetime
import logging
import os


@define(kw_only=True)
class RFDiffusionJob(BatchFoldJob):
    """Define RFDesign Hallucinate Job"""

    input_s3_uri: str
    output_s3_uri: str
    params: dict = {}
    weights_dir: str = "/database"
    output_dir: str = "outputs"
    job_name: str = field(default="SuperfoldJob" + datetime.now().strftime("%Y%m%d%s"))
    job_definition_name: str = field(default="SuperfoldJobDefinition")
    cpu: int = field(default=4)
    memory: int = field(default=16)
    gpu: int = field(default=1)

    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""
        input_filename = os.path.basename(self.input_s3_uri)

        command_list = [f"-i {self.input_s3_uri}:/opt/superfold/inputs/"]
        command_list.extend([f"-o outputs/:{self.output_s3_uri}"])
        command_list.extend(
            [
                "python run_superfold.py",
                f"/opt/superfold/inputs/{input_filename}",
                f"--alphafold_weights={self.weights_dir}",
            ]
        )
        command_list.extend([f"{key}=\'{value}\'" for key, value in self.params.items()])
        logging.info(f"Command is \n{command_list}")
        self.define_container_overrides(command_list, self.cpu, self.memory, self.gpu)
        return None
