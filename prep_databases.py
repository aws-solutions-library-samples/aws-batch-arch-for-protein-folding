from __future__ import print_function

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
# 

import boto3
from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.download_job import DownloadJob
import logging
from datetime import datetime
import urllib3

LOGGER = logging.getLogger()
LOGGER.setLevel(logging.INFO)
http = urllib3.PoolManager()

def main():
    """ Download all data to file system """

    boto_session = boto3.session.Session()
    batch_environment = BatchFoldEnvironment(boto_session = boto_session)
    job_queue_name="CPUOnDemandJobQueue"
    response = []

    download_scripts = [
        "download_alphafold_params",
        "download_bfd",
        "download_esmfold_params",
        "download_mgnify",
        "download_omegafold_params",
        "download_openfold_params",
        "download_pdb_mmcif",
        "download_pdb_seqres",
        "download_pdb70",
        "download_rfdesign_params",
        "download_small_bfd",
        "download_test",
        "download_uniprot",
        "download_uniref30",
        "download_uniref90",
        "download_diffdock_params",
        "download_rfdiffusion_params"
    ]

    for script in download_scripts:
        response.append(
            batch_environment.submit_job(
                DownloadJob(job_name=script + datetime.now().strftime("%Y%m%dT%H%M%S"), script=f"./scripts/{script}.sh"),
                job_queue_name=job_queue_name,
            )
        )

    return(response)

if __name__ == "__main__":
    main()