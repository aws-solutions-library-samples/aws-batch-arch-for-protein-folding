# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from attrs import define, field
from batchfold.batchfold_job import BatchFoldJob
from datetime import datetime
import logging
import uuid

@define
class OpenFoldJob(BatchFoldJob):
    """Define OpenFold Job"""

    bfd_database_path: str = (
        "bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    )
    config_preset: str = "model_1_ptm"
    openfold_checkpoint_path: str = "openfold_params/finetuning_ptm_2.pt"
    data_dir: str = "/database"
    data_random_seed: str = ""
    fasta_s3_uri: str = ""
    hhblits_binary_path: str = "/usr/bin/hhblits"
    hhsearch_binary_path: str = "/usr/bin/hhsearch"
    jackhmmer_binary_path: str = "/usr/bin/jackhmmer"
    kalign_binary_path: str = "/usr/bin/kalign"
    jax_param_path: str = ""
    max_template_date: str = datetime.now().strftime("%Y-%m-%d")
    mgnify_database_path: str = "mgnify/mgy_clusters_2018_12.fa"
    model_device: str = "cuda:0"
    msa_s3_uri: str = ""
    multimer_ri_gap: int = 200
    obsolete_pdbs_path: str = "pdb_mmcif/obsolete.dat"
    output_dir: str = uuid.uuid4().hex
    output_postfix: str = ""
    output_s3_uri: str = ""
    pdb70_database_path: str = "pdb70/pdb70"
    preset: str = "full_dbs"
    release_dates_path: str = ""
    save_outputs: bool = True
    skip_relaxation: bool = False
    small_bfd_database_path: str = "small_bfd/bfd-first_non_consensus_sequences.fasta"
    subtract_plddt: bool = False
    target_id: str = datetime.now().strftime("%Y%m%d%s")
    template_mmcif_dir: str = "pdb_mmcif/mmcif_files"
    trace_model: bool = False
    uniclust30_database_path: str = "uniclust30/uniclust30_2018_08/uniclust30_2018_08"
    uniref90_database_path: str = "uniref90/uniref90.fasta"
    use_precomputed_msas: bool = True
    long_sequence_inference: bool = False
    job_name: str = field(default="OpenFoldJob" + datetime.now().strftime("%Y%m%d%s"))
    job_definition_name: str = field(default="OpenFoldJobDefinition")
    cpu: int = field(default=4)
    memory: int = field(default=15)
    gpu: int = field(default=1)

    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""

        command_list = [f"-i {self.fasta_s3_uri}:{self.output_dir}/fasta/"]
        if self.use_precomputed_msas:
            command_list.extend([f"-i {self.msa_s3_uri}/jackhmmer/:{self.output_dir}/msas/{self.target_id}/"])
        command_list.extend([f"-o {self.output_dir}/predictions/:{self.output_s3_uri}"])
        command_list.extend([f"-o {self.output_dir}/timings.json:{self.output_s3_uri}/timings.json"])
        command_list.extend([
            "python3",
            "/opt/openfold/run_pretrained_openfold.py",
            f"{self.output_dir}/fasta",
            f"{self.data_dir}/{self.template_mmcif_dir}",
            f"--model_device={self.model_device}",
            f"--config_preset={self.config_preset}",
            f"--multimer_ri_gap={self.multimer_ri_gap}",
            f"--output_dir {self.output_dir}",
            f"--max_template_date={self.max_template_date}",
            f"--jackhmmer_binary_path={self.jackhmmer_binary_path}",
            f"--hhblits_binary_path={self.hhblits_binary_path}",
            f"--hhsearch_binary_path={self.hhsearch_binary_path}",
            f"--kalign_binary_path={self.kalign_binary_path}",
            f"--obsolete_pdbs_path={self.data_dir}/{self.obsolete_pdbs_path}",
        ])
        if self.use_precomputed_msas is False:
            command_list.extend(
                [
                    f"--mgnify_database_path {self.data_dir}/{self.mgnify_database_path}",
                    f"--pdb70_database_path {self.data_dir}/{self.pdb70_database_path}",
                    f"--uniclust30_database_path {self.data_dir}/{self.uniclust30_database_path}",
                    f"--uniref90_database_path {self.data_dir}/{self.uniref90_database_path}",
                ]
            )
            if self.preset == "full_dbs":
                command_list.extend(
                    [
                        f"--bfd_database_path={self.data_dir}/{self.bfd_database_path}",
                        f"--uniclust30_database_path={self.data_dir}/{self.uniclust30_database_path}",
                    ]
                )
            elif self.preset == "reduced_dbs":
                command_list.extend(
                    [
                        f"--small_bfd_database_path={self.data_dir}/{self.small_bfd_database_path}"
                    ]
                )
            else:
                raise ValueError(
                    "preset value must be either 'full_dbs' or 'reduced_dbs'"
                )
        else:
            command_list.extend(
                [f"--use_precomputed_alignments {self.output_dir}/msas/"]
            )
        if self.release_dates_path != "":
            command_list.extend([f"--release_dates_path={self.release_dates_path}"])
        if self.output_postfix != "":
            command_list.extend([f"--output_postfix={self.output_postfix}"])
        if self.data_random_seed != "":
            command_list.extend([f"--data_random_seed={self.data_random_seed}"])
        if self.openfold_checkpoint_path != "":
            command_list.extend(
                [
                    f"--openfold_checkpoint_path={self.data_dir}/{self.openfold_checkpoint_path}"
                ]
            )
        elif self.jax_param_path != "":
            command_list.extend(
                [f"--jax_param_path={self.data_dir}/{self.jax_param_path}"]
            )
        else:
            raise ValueError(
                "Please provide a value for either openfold_checkpoint_path or jax_param_path"
            )
        if self.skip_relaxation:
            command_list.extend(["--skip_relaxation"])
        if self.trace_model:
            command_list.extend(["--trace_model"])
        if self.subtract_plddt:
            command_list.extend(["--subtract_plddt"])
        if self.save_outputs:
            command_list.extend(["--save_outputs"])
        if self.long_sequence_inference:
            command_list.extend(["--long_sequence_inference"])            
        logging.info(f"Command is \n{command_list}")
        self.define_container_overrides(command_list, self.cpu, self.memory, self.gpu)
        return None