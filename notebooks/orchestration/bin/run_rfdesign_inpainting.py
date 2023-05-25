# Import required Python packages
import boto3
from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.rfdesign_job import RFDesignInpaintJob
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
    job_queue_name = "G4dnJobQueue"

    inpainting_job_name = "RFDesignInpaintingJob" + datetime.now().strftime("%Y%m%d%s")
    job_queue_name = "G4dnJobQueue"
    params = {
        "contigs": "25-35,B63-82,15-25,B119-140,0-15",
        "len": "80-115",
        "num_designs": total_num,
        "dump_all": True,
    }
    new_job = RFDesignInpaintJob(
        boto_session=boto_session,
        job_name=inpainting_job_name,
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
