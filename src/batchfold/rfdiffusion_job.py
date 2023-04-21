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
    weights_dir: str = "/database/rfdiffusion_params"
    output_dir: str = "outputs"
    job_name: str = field(
        default="RFDiffusionJob" + datetime.now().strftime("%Y%m%d%s")
    )
    job_definition_name: str = field(default="RFDiffusionJobDefinition")
    cpu: int = field(default=4)
    memory: int = field(default=16)
    gpu: int = field(default=1)

    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""

        if "inference.model_directory_path" in self.params.keys():
            raise AttributeError("Please provide the model directory path as a separate argument.") 
        elif "inference.output_prefix" in self.params.keys():
            raise AttributeError("Please provide the output prefix as a separate argument.") 
        elif "inference.input_pdb" in self.params.keys():
            raise AttributeError("Please provide the input_pdb file via the input_s3_uri parameter.") 
        input_filename = os.path.basename(self.input_s3_uri)

        command_list = [f"-i {self.input_s3_uri}:/app/RFdiffusion/inputs/"]
        command_list.extend([f"-o outputs/:{self.output_s3_uri}"])
        command_list.extend(
            [
                "python3.9 scripts/run_inference.py",
                f"inference.input_pdb=/app/RFdiffusion/inputs/{input_filename}",
                f"inference.model_directory_path={self.weights_dir}",
                f"inference.output_prefix={self.output_dir}/output",
            ]
        )
        command_list.extend([f"{key}=\'{value}\'" for key, value in self.params.items()])
        logging.info(f"Command is \n{command_list}")
        self.define_container_overrides(command_list, self.cpu, self.memory, self.gpu)
        return None
