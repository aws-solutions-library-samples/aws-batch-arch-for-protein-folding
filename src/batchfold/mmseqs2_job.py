from attrs import define
from batchfold.batchfold_job import BatchFoldJob
from datetime import datetime
import logging

@define
class MMseqs2Job(BatchFoldJob):
    """ Define MMSeqs2 MSA Creation Job """
    job_definition_name: str = "MSAJobDefinition"
    target_id: str = datetime.now().strftime("%Y%m%d%s")
    fasta_s3_uri: str = ""
    output_s3_uri: str = ""
    output_dir: str = "/tmp/msa"
    data_dir: str = "/database"
    mmseqs_database_path: str = "mmseqs_dbs"
    mmseqs_binary_path: str = "/usr/bin/mmseqs"
    hhsearch_binary_path:str = "/usr/bin/hsearch"
    uniref_db: str = "uniref30_2103_db"
    env_db: str = "colabfold_envdb_202108_db"
    pdb70_database_path: str = "pdb70/pdb70"
    
    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""

        download_string = f"aws s3 cp {self.fasta_s3_uri} {self.output_dir}/fasta/"    

        command_list = [
            f"python3 /opt/msa/create_alignments_mmseqs.py",
            f"{self.output_dir}/fasta/{self.target_id}.fasta", 
            f"{self.data_dir}/{self.mmseqs_database_path}",
            self.uniref_db,
            self.output_dir,
            self.mmseqs_binary_path,
            f"--hhsearch_binary_path {self.hhsearch_binary_path}",
        ]
        
        command_list.extend([f"--env_db {self.env_db}"]) if self.env_db else None
        command_list.extend([f"--pdb70 {self.data_dir}/{self.pdb70_database_path}"]) if self.pdb70_database_path else None
        
        upload_string = f"aws s3 cp --recursive {self.output_dir} {self.output_s3_uri}"

        command_string = download_string + " && " + " ".join(command_list) + " && " + upload_string
        logging.info(f"Command is \n{command_string}")
        self.container_overrides["command"] = [command_string]

        return None