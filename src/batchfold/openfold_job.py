from attrs import define
from batchfold.batchfold_job import BatchFoldJob
from datetime import datetime
import logging

@define
class OpenFoldJob(BatchFoldJob):
    """ Define OpenFold Job """
    target_id: str = datetime.now().strftime("%Y%m%d%s")
    template_mmcif_dir: str = "/database/pdb_mmcif/mmcif_files" 
    fasta_s3_uri: str = ""
    output_s3_uri: str = ""    
    use_precomputed_alignments: str = ""
    model_device: str = "cuda:0"
    config_preset: str = "model_1"
    jax_param_path: str = ""
    openfold_checkpoint_path: str = ""
    save_outputs: bool = False
    preset: str = "full_dbs"
    output_postfix: str = ""
    data_random_seed: str = ""
    skip_relaxation: bool = False
    multimer_ri_gap: int = 200
    output_dir: str = "data/output"
    bfd_database_path: str = "/database/bfd_database_path/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    mgnify_database_path: str = "/database/mgnify_database_path/mgy_clusters_2018_12.fa"
    pdb70_database_path: str = "/database/pdb70_database_path/pdb70"
    obsolete_pdbs_path: str = "/database/obsolete_pdbs_path/obsolete.dat"
    pdb_seqres_database_path: str = "/database/pdb_seqres_database_path/pdb_seqres.txt"
    small_bfd_database_path: str = "/database/small_bfd_database_path/bfd-first_non_consensus_sequences.fasta"
    uniclust30_database_path: str = "/database/uniclust30_database_path/uniclust30_2018_08/uniclust30_2018_08"
    uniprot_database_path: str = "/database/uniprot_database_path/uniprot.fasta"
    uniref90_database_path: str = "/database/uniref90_database_path/uniref90.fasta"
    max_template_date=datetime.now().strftime("%Y-%m-%d")

    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""

        download_string = f"aws s3 cp {self.fasta_s3_uri} data/fasta/"    
        if self.use_precomputed_alignments != "":
            download_string += f" && aws s3 cp --recursive {self.use_precomputed_alignments} data/msas/{self.target_id}/" 

        command_list = [
            f"python3 /opt/openfold/run_pretrained_openfold.py data/fasta {self.template_mmcif_dir}",
            f"--model_device={self.model_device}",
            f"--config_preset={self.config_preset}",
            f"--multimer_ri_gap={self.multimer_ri_gap}",
            f"--output_dir {self.output_dir}",
            "--save_outputs"
        ]

        if self.use_precomputed_alignments is None:
            command_list.extend([
                f"--mgnify_database_path {self.mgnify_database_path}",
                f"--pdb70_database_path {self.pdb70_database_path}",             
                f"--uniclust30_database_path /{self.uniclust30_database_path}",
                f"--uniref90_database_path {self.uniref90_database_path}"
            ])
            if self.preset == "full_dbs":
                command_list.extend([
                    f"--bfd_database_path={self.bfd_database_path}",
                    f"--uniclust30_database_path={self.uniclust30_database_path}"                       
                    ])
            elif self.preset == "reduced_dbs":
                command_list.extend([
                    f"--small_bfd_database_path={self.small_bfd_database_path}"
                ])   
            else:
                raise ValueError("preset value must be either 'full_dbs' or 'reduced_dbs'")
        else:
            command_list.extend(["--use_precomputed_alignments data/msas"])
            
        if self.skip_relaxation:
            command_list.extend(["--skip_relaxation"])

        if self.data_random_seed != "":
            command_list.extend([f"--data_random_seed={self.data_random_seed}"])

        if self.output_postfix != "":
            command_list.extend([f"--output_postfix={self.output_postfix}"])
            
        if self.openfold_checkpoint_path != "":
            command_list.extend([f"--openfold_checkpoint_path={self.openfold_checkpoint_path}"])
        elif self.jax_param_path != "":
            command_list.extend([f"--jax_param_path={self.jax_param_path}"])
        else:
            raise ValueError("Please provide a value for either openfold_checkpoint_path or jax_param_path")
        
        upload_string = f"aws s3 cp --recursive data/output/ {self.output_s3_uri}"

        command_string = download_string + " && " + " ".join(command_list) + " && " + upload_string
        logging.info(f"OpenFold command is \n{command_string}")
        self.container_overrides["command"] = [command_string]

        return None