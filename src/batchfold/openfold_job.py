from attrs import define
from batchfold.batchfold_job import BatchFoldJob
from datetime import datetime
import logging

@define
class OpenFoldJob(BatchFoldJob):
    """ Define OpenFold Job """
    job_definition_name: str = "OpenFoldJobDefinition"
    target_id: str = datetime.now().strftime("%Y%m%d%s")
    fasta_s3_uri: str = ""
    msa_s3_uri: str = ""
    output_s3_uri: str = ""
    use_precomputed_msas: bool = True
    output_dir: str = "/tmp/openfold"
    data_dir: str = "/database"
    template_mmcif_dir: str = "pdb_mmcif/mmcif_files" 
    bfd_database_path: str = "bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    mgnify_database_path: str = "mgnify/mgy_clusters_2018_12.fa"
    pdb70_database_path: str = "pdb70/pdb70"
    obsolete_pdbs_path: str = "pdb_mmcif/obsolete.dat"
    pdb_seqres_database_path: str = "pdb_seqres/pdb_seqres.txt"
    small_bfd_database_path: str = "small_bfd/bfd-first_non_consensus_sequences.fasta"
    uniclust30_database_path: str = "uniclust30/uniclust30_2018_08/uniclust30_2018_08"
    uniprot_database_path: str = "uniprot/uniprot.fasta"
    uniref90_database_path: str = "uniref90/uniref90.fasta"    
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
    
    max_template_date=datetime.now().strftime("%Y-%m-%d")

    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""

        download_string = f"aws s3 cp {self.fasta_s3_uri} {self.output_dir}/fasta/"    
        if self.use_precomputed_msas:
            download_string += f" && aws s3 cp --recursive {self.msa_s3_uri} {self.output_dir}/msas/{self.target_id}/" 

        command_list = [
            f"python3 /opt/openfold/run_pretrained_openfold.py {self.output_dir}/fasta {self.data_dir}/{self.template_mmcif_dir}",
            f"--model_device={self.model_device}",
            f"--config_preset={self.config_preset}",
            f"--multimer_ri_gap={self.multimer_ri_gap}",
            f"--output_dir {self.output_dir}/output",
            "--save_outputs"
        ]

        if self.use_precomputed_msas is False:
            command_list.extend([
                f"--mgnify_database_path {self.data_dir}/{self.mgnify_database_path}",
                f"--pdb70_database_path {self.data_dir}/{self.pdb70_database_path}",             
                f"--uniclust30_database_path {self.data_dir}//{self.uniclust30_database_path}",
                f"--uniref90_database_path {self.data_dir}/{self.uniref90_database_path}"
            ])
            if self.preset == "full_dbs":
                command_list.extend([
                    f"--bfd_database_path={self.data_dir}/{self.bfd_database_path}",
                    f"--uniclust30_database_path={self.data_dir}/{self.uniclust30_database_path}"                       
                    ])
            elif self.preset == "reduced_dbs":
                command_list.extend([
                    f"--small_bfd_database_path={self.data_dir}/{self.small_bfd_database_path}"
                ])   
            else:
                raise ValueError("preset value must be either 'full_dbs' or 'reduced_dbs'")
        else:
            command_list.extend([f"--use_precomputed_alignments {self.output_dir}/msas"])
            
        if self.skip_relaxation:
            command_list.extend(["--skip_relaxation"])

        if self.data_random_seed != "":
            command_list.extend([f"--data_random_seed={self.data_random_seed}"])

        if self.output_postfix != "":
            command_list.extend([f"--output_postfix={self.output_postfix}"])
            
        if self.openfold_checkpoint_path != "":
            command_list.extend([f"--openfold_checkpoint_path={self.data_dir}/{self.openfold_checkpoint_path}"])
        elif self.jax_param_path != "":
            command_list.extend([f"--jax_param_path={self.data_dir}/{self.jax_param_path}"])
        else:
            raise ValueError("Please provide a value for either openfold_checkpoint_path or jax_param_path")
        
        upload_string = f"aws s3 cp --recursive {self.output_dir} {self.output_s3_uri}"

        command_string = download_string + " && " + " ".join(command_list) + " && " + upload_string
        logging.info(f"Command is \n{command_string}")
        self.container_overrides["command"] = [command_string]

        return None