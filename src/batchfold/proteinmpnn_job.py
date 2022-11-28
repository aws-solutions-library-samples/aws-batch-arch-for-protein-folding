# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from attrs import define, field
from batchfold.batchfold_job import BatchFoldJob
from datetime import datetime
import logging

@define
class ProteinMPNNJob(BatchFoldJob):
    """Define ESMFold Job"""

    pdb_s3_uri: str = ""
    jsonl_s3_uri: str = ""
    output_s3_uri: str = ""
    job_definition_name: str = "ProteinMPNNJobDefinition"
    target_id: str = datetime.now().strftime("%Y%m%d%s")
    job_name: str = field(default="ProteinMPNNJob" + datetime.now().strftime("%Y%m%d%s"))
    job_definition_name: str = field(default="ProteinMPNNJobDefinition")
    cpu: int = field(default=4)
    memory: int = field(default=15)
    gpu: int = field(default=0)

    ca_only: bool = field(default=False)
    path_to_model_weights: str = field(default="")
    model_name: str = field(default="v_48_020")
    seed: int = field(default=0)
    save_score: int = field(default=0)
    save_probs: int = field(default=0)
    score_only: int = field(default=0)
    conditional_probs_only: int = field(default=0)
    conditional_probs_only_backbone: int = field(default=0)
    unconditional_probs_only: int = field(default=0)
    backbone_noise: float = field(default=0.00)
    num_seq_per_target: int = field(default=1)
    batch_size: int = field(default=1)
    max_length: int = field(default=200000)
    sampling_temp: str = field(default="0.1")
    pdb_path_chains: str = field(default="")
    chain_id_jsonl: str = field(default="")
    fixed_positions_jsonl: str = field(default="")
    omit_AAs: str = field(default="X")
    bias_AA_jsonl: str = field(default="")
    bias_by_res_jsonl: str = field(default="")
    omit_AA_jsonl: str = field(default="")
    pssm_jsonl: str = field(default="")
    pssm_multi: float = field(default=0.00)
    pssm_threshold: float = field(default=0.00)
    pssm_log_odds_flag: int = field(default=0)
    pssm_bias_flag: int = field(default=0)
    tied_positions_jsonl: str = field(default="")

    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""

        if self.pdb_s3_uri:
            input_command = f"-i {self.pdb_s3_uri}:input/{self.target_id}.pdb"
        elif self.jsonl_s3_uri:
            input_command = f"-i {self.jsonl_s3_uri}:input/{self.target_id}.jsonl"
        else:
            raise ValueError(
                "Please specify a value for either pdb_s3_uri or pdb_s3_uri"
            )

        command_list = [input_command]
        command_list.extend([f"-o output/:{self.output_s3_uri}"])
        command_list.extend(["python protein_mpnn_run.py"])
        command_list.extend([f"--out_folder=output"])
        
        if self.pdb_s3_uri:
            command_list.extend([f"--pdb_path=input/{self.target_id}.pdb"])
            command_list.extend([f"--pdb_path_chains='{self.pdb_path_chains}'"])

        elif self.jsonl_s3_uri:
            command_list.extend([f"--jsonl_s3_uri=input/{self.target_id}.jsonl"])

        if self.ca_only:
            command_list.extend([f"--ca_only={self.ca_only}"])
        if self.path_to_model_weights:
            command_list.extend([f"--path_to_model_weights={self.path_to_model_weights}"])
        if self.model_name:
            command_list.extend([f"--model_name={self.model_name}"])
        if self.seed:
            command_list.extend([f"--seed={self.seed}"])
        if self.save_score:
            command_list.extend([f"--save_score={self.save_score}"])
        if self.save_probs:
            command_list.extend([f"--save_probs={self.save_probs}"])
        if self.score_only:
            command_list.extend([f"--score_only={self.score_only}"])
        if self.conditional_probs_only:
            command_list.extend([f"--conditional_probs_only={self.conditional_probs_only}"])
        if self.conditional_probs_only_backbone:
            command_list.extend([f"--conditional_probs_only_backbone={self.conditional_probs_only_backbone}"])
        if self.unconditional_probs_only:
            command_list.extend([f"--unconditional_probs_only={self.unconditional_probs_only}"])
        if self.backbone_noise:
            command_list.extend([f"--backbone_noise={self.backbone_noise}"])
        if self.num_seq_per_target:
            command_list.extend([f"--num_seq_per_target={self.num_seq_per_target}"])
        if self.batch_size:
            command_list.extend([f"--batch_size={self.batch_size}"])
        if self.max_length:
            command_list.extend([f"--max_length={self.max_length}"])
        if self.sampling_temp:
            command_list.extend([f"--sampling_temp={self.sampling_temp}"])
        if self.chain_id_jsonl:
            command_list.extend([f"--chain_id_jsonl={self.chain_id_jsonl}"])
        if self.fixed_positions_jsonl:
            command_list.extend([f"--fixed_positions_jsonl={self.fixed_positions_jsonl}"])
        if self.omit_AAs:
            command_list.extend([f"--omit_AAs={self.omit_AAs}"])
        if self.bias_AA_jsonl:
            command_list.extend([f"--bias_AA_jsonl={self.bias_AA_jsonl}"])
        if self.bias_by_res_jsonl:
            command_list.extend([f"--bias_by_res_jsonl={self.bias_by_res_jsonl}"])
        if self.omit_AA_jsonl:
            command_list.extend([f"--omit_AA_jsonl={self.omit_AA_jsonl}"])
        if self.pssm_jsonl:
            command_list.extend([f"--pssm_jsonl={self.pssm_jsonl}"])
        if self.pssm_multi:
            command_list.extend([f"--pssm_multi={self.pssm_multi}"])
        if self.pssm_threshold:
            command_list.extend([f"--pssm_threshold={self.pssm_threshold}"])
        if self.pssm_log_odds_flag:
            command_list.extend([f"--pssm_log_odds_flag={self.pssm_log_odds_flag}"])
        if self.pssm_bias_flag:
            command_list.extend([f"--pssm_bias_flag={self.pssm_bias_flag}"])
        if self.tied_positions_jsonl:
            command_list.extend([f"--tied_positions_jsonl={self.tied_positions_jsonl}"])      
        
        logging.info(f"Command is \n{command_list}")
        self.define_container_overrides(command_list, self.cpu, self.memory, self.gpu)
        return None