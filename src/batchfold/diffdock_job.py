# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from attrs import define, field
from batchfold.batchfold_job import BatchFoldJob
from datetime import datetime
import logging

@define
class DiffDockJob(BatchFoldJob):
    """Define DiffDock Job"""

    protein_s3_uri: str = ""
    ligand_s3_uri: str = ""
    output_s3_uri: str = ""
    complex_name: str = datetime.now().strftime("%Y%m%d%s")
    protein_path: str = "input/protein.pdb"
    ligand_description: str = "input/ligand.sdf"
    out_dir: str = "output"
    save_visualisation: bool = True
    samples_per_complex: int = 10
    model_dir: str = "workdir/paper_score_model"
    ckpt: str = "best_ema_inference_epoch_model.pt"
    confidence_model_dir: str = "workdir/paper_confidence_model"
    confidence_ckpt: str = "best_model_epoch75.pt"
    batch_size: int = 32
    no_final_step_noise: bool = False
    inference_steps: int = 20
    actual_steps: int = 20
    so_distribution_cache_dir: str = "/database/diffdock_params"
    esm_param_model_dir: str = "/database/diffdock_params"
    job_name: str = field(default="DiffDockJob" + datetime.now().strftime("%Y%m%d%s"))
    job_definition_name: str = field(default="DiffDockJobDefinition")
    cpu: int = field(default=4)
    memory: int = field(default=15)
    gpu: int = field(default=1)

    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""
        command_list = [f"-i {self.protein_s3_uri}:{self.protein_path}"]
        if self.ligand_description.endswith(".sdf"):
            command_list.extend([f"-i {self.ligand_s3_uri}:{self.ligand_description}"])
        command_list.extend([f"-o {self.out_dir}/:{self.output_s3_uri}"])
        command_list.extend([
            "python inference.py",
            f"--complex_name={self.complex_name}",
            f"--protein_path={self.protein_path}",
            f"--ligand_description='{self.ligand_description}'",
            f"--out_dir={self.out_dir}",
            f"--samples_per_complex={self.samples_per_complex}",
            f"--model_dir={self.model_dir}",
            f"--confidence_model_dir={self.confidence_model_dir}",
            f"--confidence_ckpt={self.confidence_ckpt}",
            f"--batch_size={self.batch_size}",
            f"--inference_steps={self.inference_steps}",
            f"--actual_steps={self.actual_steps}",
            f"--so_distribution_cache_dir={self.so_distribution_cache_dir}",
            f"--esm_param_model_dir={self.esm_param_model_dir}",
        ])

        if self.save_visualisation:
            command_list.extend([f"--save_visualisation"])
        if self.no_final_step_noise:
            command_list.extend(["--no_final_step_noise"])

        logging.info(f"Command is \n{command_list}")
        self.define_container_overrides(command_list, self.cpu, self.memory, self.gpu)
        return None