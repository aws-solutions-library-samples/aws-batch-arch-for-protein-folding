from attrs import define
from batchfold.batchfold_job import BatchFoldJob
from datetime import datetime
import logging

@define
class AlphaFold2Job(BatchFoldJob):
    """ Define AlphaFold 2 Job """
    target_id: str = datetime.now().strftime("%Y%m%d%s")
    fasta_s3_uri: str = ""
    msa_s3_uri: str = ""
    output_s3_uri: str = ""
    data_dir: str = "/data"
    output_dir: str = "/data/output"
    bfd_database_path: str = "/database/bfd_database_path/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    mgnify_database_path: str = "/database/mgnify_database_path/mgy_clusters_2018_12.fa"
    pdb70_database_path: str = "/database/pdb70_database_path/pdb70"
    obsolete_pdbs_path: str = "/database/obsolete_pdbs_path/obsolete.dat"
    template_mmcif_dir: str = "/database/template_mmcif_dir/mmcif_files"
    pdb_seqres_database_path: str = "/database/pdb_seqres_database_path/pdb_seqres.txt"
    small_bfd_database_path: str = "/database/small_bfd_database_path/bfd-first_non_consensus_sequences.fasta"
    uniclust30_database_path: str = "/database/uniclust30_database_path/uniclust30_2018_08/uniclust30_2018_08"
    uniprot_database_path: str = "/database/uniprot_database_path/uniprot.fasta"
    uniref90_database_path: str = "/database/uniref90_database_path/uniref90.fasta"
    max_template_date: str = datetime.now().strftime("%Y-%m-%d")
    db_preset: str = "reduced_dbs"
    model_preset: str = "monomer"
    benchmark: bool = False
    use_precomputed_msa: bool = False
    features_paths: str = ""
    run_features_only: bool = False
    logtostderr: str = True
    run_relax: bool = True
    use_gpu_relax: bool = True
    num_multimer_predictions_per_model: int = 1   

    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""

        download_string = f"aws s3 cp {self.fasta_s3_uri} data/fasta/"    
        if self.use_precomputed_msa:
            download_string += f" && aws s3 cp --recursive {self.msa_s3_uri} data/{self.target_id}/msas/" 

        command_list = [
            f"/app/run_alphafold.sh",
            f"--fasta_paths=/data/fasta/{self.target_id}.fasta",
            "--data_dir=/database",
            "--output_dir=/data",
            "--template_mmcif_dir=/mnt/pdb_mmcif/mmcif_files",
            "--obsolete_pdbs_path=/mnt/pdb_mmcif/obsolete.dat",
            # "--uniref90_database_path=/database/uniref90/uniref90.fasta",
            # "--mgnify_database_path=/mnt/mgnify/mgy_clusters_2018_12.fa",            
            # "--small_bfd_database_path=/mnt/small_bfd/bfd-first_non_consensus_sequences.fasta",
            # "--pdb70_database_path=/mnt/pdb70/pdb70",
            # f"--db_preset={self.db_preset}",
            f"--max_template_date={self.max_template_date}",
            f"--model_preset={self.model_preset}",
            f"--use_precomputed_msas={self.use_precomputed_msa}",
            f"--benchmark={self.benchmark}",
            f"--run_relax={self.run_relax}",
            f"--use_gpu_relax={self.use_gpu_relax}",
            f"--num_multimer_predictions_per_model={self.num_multimer_predictions_per_model}"
        ]

        if self.use_precomputed_msa is None:
            command_list.extend([      
                f"--db_preset={self.db_preset}",          
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
        
        upload_string = f"aws s3 cp --recursive data/output/ {self.output_s3_uri}"

        command_string = download_string + " && " + " ".join(command_list) + " && " + upload_string
        logging.info(f"AlphaFold2 command is \n{command_string}")
        self.container_overrides["command"] = [command_string]

        return None