from attrs import define
from batchfold.batchfold_job import BatchFoldJob
from datetime import datetime
from pathlib import Path
import logging

@define
class DownloadJob(BatchFoldJob):
    """ Define Jackhmmer MSA Creation Job """
    job_definition_name: str = "MSAJobDefinition"
    script: str = "./scripts/download_test.sh"
    data_dir: str = "/database"
    
    def __attrs_post_init__(self) -> None:
        """Override default BatchFoldJob command"""

        self.job_name = Path(self.script).stem + "-" + datetime.now().strftime("%Y%m%dT%H%M%S"),

        command_string = f"{self.script} {self.data_dir}"
        
        logging.info(f"Command is \n{command_string}")
        self.container_overrides["command"] = [command_string]

        return None