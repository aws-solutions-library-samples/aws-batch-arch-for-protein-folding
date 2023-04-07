# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from attrs import define, field
from batchfold.batchfold_job import BatchFoldJob
from datetime import datetime
import logging
import uuid

@define
class ESMFoldJob(BatchFoldJob):
    """Define ESMFold Job"""

    fasta_s3_uri: str = ""
    output_s3_uri: str = ""
    target_id: str = datetime.now().strftime("%Y%m%d%s")
    num_recycles: int = 4
    max_tokens_per_batch: int = 1024
    chunk_size: int = None
    cpu_only: bool = False
    cpu_offload: bool = False
    model_dir: str = "/database/esmfold_params/hub"
    job_name: str = field(default="ESMFoldJob" + datetime.now().strftime("%Y%m%d%s"))
    job_definition_name: str = field(default="ESMFoldJobDefinition")
    cpu: int = field(default=8)
    memory: int = field(default=31)
    gpu: int = field(default=1)

    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""
        command_list = [f"-i {self.fasta_s3_uri}:input/{self.target_id}.fasta"]
        command_list.extend([f"-o output/:{self.output_s3_uri}"])
        command_list.extend([
            "python scripts/esmfold_inference.py",
            f"--fasta=input/{self.target_id}.fasta",
            f"--pdb=output",
            f"--num-recycles={self.num_recycles}",
            f"--max-tokens-per-batch={self.max_tokens_per_batch}",
            f"--model-dir={self.model_dir}",
        ])

        if self.chunk_size:
            command_list.extend([f"--chunk-size={self.chunk_size}"])
        if self.cpu_only:
            command_list.extend(["--cpu_only"])
        if self.cpu_offload:
            command_list.extend(["--cpu_offload"])

        logging.info(f"Command is \n{command_list}")
        self.define_container_overrides(command_list, self.cpu, self.memory, self.gpu)
        return None