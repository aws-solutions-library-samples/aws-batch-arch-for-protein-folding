# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from attrs import define
from batchfold.batchfold_job import BatchFoldJob
from datetime import datetime
import logging


@define
class JackhmmerJob(BatchFoldJob):
    """Define Jackhmmer MSA Creation Job"""

    fasta_s3_uri: str = ""
    output_s3_uri: str = ""
    output_dir: str = "/tmp/msa"
    data_dir: str = "/database"
    uniref90_database_path: str = "uniref90/uniref90.fasta"
    mgnify_database_path: str = "mgnify/mgy_clusters_2018_12.fa"
    template_mmcif_dir: str = "pdb_mmcif/mmcif_files"
    max_template_date: str = datetime.now().strftime("%Y-%m-%d")
    obsolete_pdbs_path: str = "pdb_mmcif/obsolete.dat"
    job_definition_name: str = "MSAJobDefinition"
    pdb70_database_path: str = "pdb70/pdb70"
    bfd_database_path: str = (
        "bfd_database_path/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    )
    small_bfd_database_path: str = (
        "bfd_database_path/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    )
    uniclust30_database_path: str = (
        "uniclust30_database_path/uniclust30_2018_08/uniclust30_2018_08"
    )
    uniprot_database_path: str = "uniprot_database_path/uniprot.fasta"
    pdb_seqres_database_path: str = "pdb_seqres_database_path/pdb_seqres.txt"
    db_preset: str = "full_dbs"
    model_preset: str = "monomer"
    target_id: str = datetime.now().strftime("%Y%m%d%s")

    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""

        download_string = f"aws s3 cp {self.fasta_s3_uri} {self.output_dir}/fasta/"

        command_list = [
            f"python3 /opt/msa/create_alignments.py",
            f"--fasta_paths {self.output_dir}/fasta",
            f"--output_dir {self.output_dir}/output",
            f"--cpu {self.cpu}",
            f"--uniref90_database_path {self.data_dir}/{self.uniref90_database_path}",
            f"--mgnify_database_path {self.data_dir}/{self.mgnify_database_path}",
            f"--template_mmcif_dir {self.data_dir}/{self.template_mmcif_dir}",
            f"--max_template_date {self.data_dir}/{self.template_mmcif_dir}",
            f"--obsolete_pdbs_path {self.data_dir}/{self.obsolete_pdbs_path}",
            f"--db_preset {self.db_preset}",
            f"--model_preset {self.model_preset}",
        ]

        if self.db_preset == "reduced_dbs":
            command_list.extend(
                [
                    f"--small_bfd_database_path {self.data_dir}/{self.small_bfd_database_path}"
                ]
            )
        else:
            command_list.extend(
                [
                    f"--bfd_database_path {self.data_dir}/{self.bfd_database_path}",
                    f"--uniclust30_database_path {self.data_dir}/{self.uniclust30_database_path}",
                ]
            )

        if self.model_preset == "multimer":
            command_list.extend(
                [
                    f"--pdb_seqres_database_path {self.data_dir}/{self.pdb_seqres_database_path}",
                    f"--uniprot_database_path {self.data_dir}/{self.uniprot_database_path}",
                ]
            )
        else:
            command_list.extend(
                [f"--pdb70_database_path {self.data_dir}/{self.pdb70_database_path}"]
            )

        upload_string = f"aws s3 cp --recursive {self.output_dir}/{self.target_id}/msas/ {self.output_s3_uri}/jackhmmer && aws s3 cp {self.output_dir}/{self.target_id}/features.pkl {self.output_s3_uri}/jackhmmer"
        command_string = (
            download_string + " && " + " ".join(command_list) + " && " + upload_string
        )
        logging.info(f"Command is \n{command_string}")
        self.container_overrides["command"] = [command_string]

        return None
