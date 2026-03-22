# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from attrs import define, field
from batchfold.batchfold_job import BatchFoldJob
from datetime import datetime
import logging
import uuid

@define
class MMSeqs2Job(BatchFoldJob):
    """Define MMSeqs2 MSA Creation Job"""

    fasta_s3_uri: str = ""
    output_s3_uri: str = ""
    output_dir: str = uuid.uuid4().hex
    data_dir: str = "/database"
    target_id: str = datetime.now().strftime("%Y%m%d%s")

    # Database selection
    databases: list = ["uniref90"]
    databases_paths_overwrite: dict = None

    # MMSeqs2 search parameters
    evalue: float = 0.1
    cov: float = 0.5
    cov_mode: int = 2
    max_seqs: int = 10000
    sensitivity: float = 7.5
    num_iterations: int = 3
    output_format: str = "a3m"
    write_hits_table: bool = True

    # AWS Batch resource configuration
    job_name: str = field(default="MMSeqs2Job" + datetime.now().strftime("%Y%m%d%s"))
    job_definition_name: str = field(default="MMSeqs2JobDefinition")
    cpu: int = field(default=16)
    memory: int = field(default=31)
    gpu: int = field(default=0)

    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""

        # Auto-select GPU job definition if GPU is requested
        if self.gpu > 0 and self.job_definition_name == "MMSeqs2JobDefinition":
            self.job_definition_name = "MMSeqs2GPUJobDefinition"

        # Build command list following run.sh convention
        command_list = [f"-i {self.fasta_s3_uri}:{self.output_dir}/fasta/{self.target_id}.fasta"]
        command_list.extend([f"-o {self.output_dir}/:{self.output_s3_uri}/mmseqs2/"])

        # Build database paths (dbname:dbpath format)
        db_paths = []
        for db in self.databases:
            if self.databases_paths_overwrite and db in self.databases_paths_overwrite:
                path = self.databases_paths_overwrite[db]
            else:
                path = f"{db}/{db}"
            db_paths.append(f"{db}:{path}")

        # Build the MMSeqs2 command
        command_list.extend([
            "bash", "/opt/mmseqs/run_mmseqs_msa.sh",
            f"--query_fasta {self.output_dir}/fasta/{self.target_id}.fasta",
            f"--run_dir {self.output_dir}",
            f"--databases_paths {','.join(db_paths)}",
            f"--threads {self.cpu}",
            f"--evalue {self.evalue}",
            f"--cov {self.cov}",
            f"--cov_mode {self.cov_mode}",
            f"--max_seqs {self.max_seqs}",
            f"--s {self.sensitivity}",
            f"--num_iterations {self.num_iterations}",
            f"--output_format {self.output_format}",
            f"--write_hits_table {int(self.write_hits_table)}",
        ])

        # Add GPU flag if GPU resources are requested
        if self.gpu > 0:
            command_list.append("--gpu 1")

        logging.info(f"Command is \n{command_list}")
        self.define_container_overrides(command_list, self.cpu, self.memory, self.gpu)
        return None
