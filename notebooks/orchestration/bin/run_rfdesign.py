# Import required Python packages
import boto3
from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.rfdesign_job import RFDesignHallucinateJob
from datetime import datetime
import logging
import sys
import argparse

logging.basicConfig(stream=sys.stderr, level=logging.INFO)


def main(args):

    input_s3_uri = args["input_s3_uri"]
    output_s3_uri = args["output_s3_uri"]
    num_sequences_to_generate = args["num_sequences_to_generate"]

    # Create AWS clients
    boto_session = boto3.session.Session()
    batch_environment = BatchFoldEnvironment(boto_session=boto_session)

    total_num = num_sequences_to_generate
    mask = "25-35,B63-82,15-25,B119-140,0-15"
    hallucinate_job_prefix = "RFDesignHallucinateJob" + datetime.now().strftime(
        "%Y%m%d%s"
    )
    job_queue_name = "G4dnJobQueue"

    job_name = f"{hallucinate_job_prefix}_0"
    params = {
        "mask": mask,
        "steps": "g10",
        "num": total_num,
        "start_num": 0,
        "w_rog": 1,
        "rog_thresh": 16,
        "w_rep": 2,
        "rep_pdb": "input/pdl1.pdb",
        "rep_sigma": 4,
        "save_pdb": True,
        "track_step": 10,
    }

    new_job = RFDesignHallucinateJob(
        boto_session=boto_session,
        job_name=job_name,
        target_id="4ZQK",
        input_s3_uri=input_s3_uri,
        output_s3_uri=output_s3_uri,
        pdb="input/pd1.pdb",
        params=params,
    )

    submission = batch_environment.submit_job(new_job, job_queue_name)
    print(submission.job_id)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Parse the options")
    parser.add_argument(
        "--input_s3_uri",
        dest="input_s3_uri",
        default=None,
        help="input s3 uri (default: None)",
    )
    parser.add_argument(
        "--output_s3_uri",
        dest="output_s3_uri",
        default=None,
        help="output_s3_uri (default: None)",
    )
    parser.add_argument(
        "--num_sequences_to_generate",
        dest="num_sequences_to_generate",
        default=1,
        help="number of sequences for rfdesign to generate (default: 1)",
    )

    args = parser.parse_args()
    args = vars(args)

    main(args)
