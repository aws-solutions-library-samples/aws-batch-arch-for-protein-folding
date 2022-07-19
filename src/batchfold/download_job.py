from attrs import define
from batchfold.batchfold_job import BatchFoldJob
from datetime import datetime
from typing import List, Dict
import logging

@define
class DownloadJob(BatchFoldJob):
    """ Define Jackhmmer MSA Creation Job """
    job_definition_name: str = "MSAJobDefinition"
    script: str = "./scripts/download_test.sh"
    data_dir: str = "/database"
    
    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""

        command_string = f"{self.script} {self.data_dir}"
        
        logging.info(f"Command is \n{command_string}")
        self.container_overrides["command"] = [command_string]

        return None