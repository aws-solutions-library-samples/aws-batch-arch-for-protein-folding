from attrs import define
from batchfold.batchfold_job import BatchFoldJob
from datetime import datetime
import logging

@define
class JackhmmerJob(BatchFoldJob):
    """ Define Jackhmmer MSA Creation Job """
    job_definition_name: str = "MSAJobDefinition"
    target_id: str = datetime.now().strftime("%Y%m%d%s")
    fasta_s3_uri: str = ""
    output_s3_uri: str = ""
    output_dir: str = "/tmp/msa"
    data_dir: str = "/database"
    bfd_database_path: str = "bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    mgnify_database_path: str = "mgnify/mgy_clusters_2018_12.fa"
    pdb70_database_path: str = "pdb70/pdb70"
    uniclust30_database_path: str = "uniclust30/uniclust30_2018_08/uniclust30_2018_08"
    uniprot_database_path: str = "uniprot/uniprot.fasta"
    uniref90_database_path: str = "uniref90/uniref90.fasta"    
    use_small_bfd: bool = False
    
    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""

        download_string = f"aws s3 cp {self.fasta_s3_uri} {self.output_dir}/fasta/"    

        command_list = [
            f"python3 /opt/msa/create_alignments.py {self.output_dir}/fasta",
            f"--output_dir {self.output_dir}/output",
            f"--cpu {self.cpu}",
        ]
        
        command_list.extend([f"--mgnify_database_path {self.data_dir}/{self.mgnify_database_path}"]) if self.mgnify_database_path else None
        command_list.extend([f"--pdb70_database_path {self.data_dir}/{self.pdb70_database_path}"]) if self.pdb70_database_path else None
        command_list.extend([f"--uniclust30_database_path {self.data_dir}/{self.uniclust30_database_path}"]) if self.uniclust30_database_path else None
        command_list.extend([f"--uniref90_database_path {self.data_dir}/{self.uniref90_database_path}"]) if self.uniref90_database_path else None
        command_list.extend([f"--bfd_database_path {self.data_dir}/{self.bfd_database_path}"]) if self.bfd_database_path else None
        if self.use_small_bfd is False:
            command_list.extend([f"--uniclust30_database_path {self.data_dir}/{self.uniclust30_database_path}"]) if self.uniclust30_database_path else None
        else:
            command_list.extend([
                f"--use_small_bfd {self.use_small_bfd}"                       
                ])
        
        upload_string = f"perl /tmp/hh-suite/scripts/reformat.pl a3m sto '{self.output_dir}/output/alignments/{self.target_id}/*.a3m' .sto"
        upload_string += f" && aws s3 cp --recursive {self.output_dir}/output/alignments/{self.target_id} {self.output_s3_uri}/jackhmmer"
        command_string = download_string + " && " + " ".join(command_list) + " && " + upload_string
        logging.info(f"Command is \n{command_string}")
        self.container_overrides["command"] = [command_string]

        return None